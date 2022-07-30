#lang typed/racket/base

;; SPDX-License-Identifier: LGPL-3.0-or-later
;;
;; tiling.rkt -- functions to determine the list of geoids which cover regions
;; on the Earth surface.
;;
;; This file is part of geoid -- work efficiently with geographic data
;; Copyright (c) 2022 Alex Harsányi <AlexHarsanyi@gmail.com>
;;
;; This program is free software: you can redistribute it and/or modify it
;; under the terms of the GNU Lesser General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or (at your
;; option) any later version.
;;
;; This program is distributed in the hope that it will be useful, but WITHOUT
;; ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
;; FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
;; License for more details.
;;
;; You should have received a copy of the GNU Lesser General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.


(require "vmath.rkt"
         "geoid.rkt"
         racket/match
         racket/math
         racket/set
         racket/list
         racket/bool
         racket/format
         racket/flonum
         math/flonum
         typed/racket/class)

(provide (all-defined-out))


;;.............................................................. helpers ....

;; Return the squared chord length subtended by two points A and B on the unit
;; sphere.  A and B are unit vectors and the returned value is between 0 (for
;; identical points) and 4 (opposite points on the unit sphere)
(: chlen²-from-points (-> Vector3 Vector3 Flonum))
(define (chlen²-from-points a b)
  ;; The squared distance can exceed 4.0 due to roundoff errors, but we make
  ;; sure that it is a maximum of 4.0
  (min 4.0 (v3len² (v3- a b))))

;; Return the squared chord length subtended by the R angle (in radians) on
;; the unit circle
(: radians->chlen² (-> Flonum Flonum))
(define (radians->chlen² r)
  (sqr (* 2.0 (sin (* 0.5 (min pi r))))))

;; Return the angle (in radians) subtended by the squared chord length C on
;; the unit circle.
(: chlen²->radians (-> Flonum Flonum))
(define (chlen²->radians c)
  (* 2.0 (assert (asin (* 0.5 (sqrt c))) flonum?)))

;; Return a unit vector for a latitude/longitude point.
(: ->unit-vector (-> (Vector Real Real) Vector3))
(define (->unit-vector lat-lon)
  (match-define (vector lat lon) lat-lon)
  (define-values (x y z) (lat-lng->unit-vector lat lon))
  (flvector x y z))


;;............................................................. segments ....

;; A segment defined on the unit sphere: START and END are the unit vectors at
;; the two ends of the segment and EDGE is the normal vector to the plane
;; defined by the two points on the unit sphere.
(struct segment ([start : Vector3]
                 [end : Vector3]
                 [edge : Vector3]) #:transparent)

;; Construct a segment from two points (the edge can be calculated)
(: make-segment (-> Vector3 Vector3 segment))
(define (make-segment start end)
  (let ([edge (v3unit (robust-cross-product start end))])
    (segment start end edge)))

;; Return true if `s` is a valid segment, that is, it has a non NaN edge --
;; trying to construct a segment with the same start and end point produces
;; invalid segments.
(: valid-segment? (-> segment Boolean))
(define (valid-segment? s)
  (define edge (segment-edge s))
  (not (or (nan? (flvector-ref edge 0))
           (nan? (flvector-ref edge 1))
           (nan? (flvector-ref edge 2)))))

;; Return true if the SEGMENT is entirely inside the CELL.  It is in the cell
;; if both the start and ending points are inside the cell.
(: segment-inside-cell? (-> segment cell Boolean))
(define (segment-inside-cell? segment cell)
  (and (point-inside-cell? (segment-start segment) cell)
       (point-inside-cell? (segment-end segment) cell)))

;; Return true if the segment S intersects the cell C.  We only return true if
;; the segment actually crosses the cell edge.
(: segment-intersects-cell?/strict (-> segment cell Boolean))
(define (segment-intersects-cell?/strict s c)
  (define start-inside? (point-inside-cell? (segment-start s) c))
  (define end-inside? (point-inside-cell? (segment-end s) c))
  (and
   ;; In strict mode, the start and end can't both be inside the cell.
   (not (and start-inside? end-inside?))
   (match-let ([(list c0 c1 c2 c3) (cell-corners c)]
               [(list e0 e1 e2 e3) (cell-edges/raw c)])
     (or (segment-or-vertex-crossing? s (segment c0 c1 e0))
         (segment-or-vertex-crossing? s (segment c1 c2 e1))
         (segment-or-vertex-crossing? s (segment c2 c3 e2))
         (segment-or-vertex-crossing? s (segment c3 c0 e3))))))

;; Return true if the segment S intersects the cell C.  We also return true if
;; the segment is entirely inside the cell.
(: segment-intersects-cell?/relaxed (-> segment cell Boolean))
(define (segment-intersects-cell?/relaxed s c)
  (or
   (point-inside-cell? (segment-start s) c)
   (point-inside-cell? (segment-end s) c)
   (match-let ([(list c0 c1 c2 c3) (cell-corners c)]
               [(list e0 e1 e2 e3) (cell-edges/raw c)])
     (or (segment-or-vertex-crossing? s (segment c0 c1 e0))
         (segment-or-vertex-crossing? s (segment c1 c2 e1))
         (segment-or-vertex-crossing? s (segment c2 c3 e2))
         (segment-or-vertex-crossing? s (segment c3 c0 e3))))))

;; Convert a track (a list of latitude/longitude points) into a list of
;; segments, each segment connecting two adjacent points on the track.  If
;; CLOSED? is #t, a segment which connects the last point on the track with
;; the first one is also created.
(: track->segments (->* ((Listof (Vector Real Real)))
                        (Boolean)
                        (Listof segment)))
(define (track->segments track [closed? #f])
  (if (null? track)
      '()
      (let loop : (Listof segment)
           ([p1 : Vector3 (->unit-vector (car track))]
            [remaining : (Listof (Vector Real Real)) (cdr track)]
            [segments : (Listof segment) '()])
        (if (null? remaining)
            (reverse
             (if closed?
                 (let ([s (make-segment p1 (->unit-vector (car track)))])
                   (if (valid-segment? s)
                       (cons s segments)
                       segments))
                 segments))
            (let ([p2 (->unit-vector (car remaining))])
              (loop p2
                    (cdr remaining)
                    (let ([s (make-segment p1 p2)])
                      (if (valid-segment? s)
                          (cons s segments)
                          segments))))))))

;; Return a geoid at LEVEL used to start a search on segment S
(: segment->start-geoid (-> segment Integer Integer))
(define (segment->start-geoid s level)
  (define start-point (segment-start s))
  (unit-vector->geoid (flvector-ref start-point 0)
                      (flvector-ref start-point 1)
                      (flvector-ref start-point 2)
                      level))

;; Return the list of geoids which contain the entire segment S.  All the
;; geoids will be at the specified LEVEL.
(: geoid-covering-for-segment (-> segment Integer (Listof Integer)))
(define (geoid-covering-for-segment s level)
  (define start-geoid (segment->start-geoid s level))
  (let loop : (Listof Integer)
       ([intersects : (Setof Integer) (set)]
        [remaining : (Listof Integer) (list start-geoid)]
        [considered : (Setof Integer) (set)])
    (if (null? remaining)
        (set->list intersects)
        (let ([geoid (car remaining)]
              [remaining (cdr remaining)])
          (if (or (set-member? intersects geoid)
                  (set-member? considered geoid))
              (loop intersects remaining considered)
              (let ([cell (make-cell geoid)])
                (if (segment-intersects-cell?/relaxed s cell)
                    (loop (set-add intersects geoid)
                          (append (adjacent-geoids geoid) remaining)
                          considered)
                    (loop intersects
                          remaining
                          (set-add considered geoid)))))))))

;; Returns +1 if the two segments cross each other at a point which is
;; interior to both edges, 0 if any two vertices from different edges are the
;; same and -1 otherwise.
(: segment-crossing-sign (-> segment segment Integer))
(define (segment-crossing-sign s1 s2)
  (define a (segment-start s1))
  (define b (segment-end s1))
  (define a-cross-b (segment-edge s1))
  (define c (segment-start s2))
  (define d (segment-end s2))
  (define c-cross-d (segment-edge s2))
  ;; Note: arguments to the sign/c can be rotated withtout changing the
  ;; result.
  (define acb (sign/c a c b (cross-product a c)))
  (define bda (sign/c a b d a-cross-b))
  (if (or (zero? acb) (zero? bda))
      0                       ; we would need to use more precise calculations
      (if (= acb (- bda))
          -1
          (let ([cbd (sign/c c b d (cross-product c b))])
            (if (zero? cbd)
                0
                (if (= acb (- cbd))
                    -1
                    (let ([dac (sign/c c d a c-cross-d)])
                      (if (zero? dac)
                          0
                          (if (= acb (- dac))
                              -1
                              1)))))))))

;; Return an arbitrary, but consistent, direction for the vector A.  This
;; direction is on the plane orthogonal to A, and it is used by
;; `angle-contains-vertex?` and `vertex-crossing?` to classify vertices as
;; in/out of an intersection.  The convention in this file (and the Google S2
;; library) is that a vertex A "belongs" only to the segment which is the
;; reference direction of "A" -- this helps with classifying loop vertices as
;; being inside or outside a loop.
(: reference-direction (-> Vector3 Vector3))
(define (reference-direction a)
  (define k
    (let ([k0 (- (v3-largest-abs-component a) 1)])
      (if (< k0 0) 2 k0)))
  (define temp (flvector 0.012 0.0053 0.00457))
  (flvector-set! temp k 1.0)
  (v3unit (cross-product a temp)))

;; See the AngleContainsVertex function in s2edge_crossing.h.  Return true if
;; the vertex B is contained in the angle A B C.  This is an arbitrary
;; definition, but consistent with `vertex-crossing?`: the angle A B C
;; contains vertex B if the reference direction for the vector B is inside the
;; angle when traversing it in CCW order.
(: angle-contains-vertex? (-> Vector3 Vector3 Vector3 Boolean))
(define (angle-contains-vertex? a b c)
  (not (ordered-ccw (reference-direction b) c a b)))

;; Return true if the segments S1 and S2 are "crossing" -- this function
;; assumes that S1 and S1 share a common start/end vertex and the "crossing"
;; definition is made to be consistent with `angle-contains-vertex?`
(: vertex-crossing? (-> segment segment Boolean))
(define (vertex-crossing? s1 s2)
  (define a (segment-start s1))
  (define b (segment-end s1))
  (define c (segment-start s2))
  (define d (segment-end s2))
  (define result
    (cond
      ((or (equal? a b) (equal? c d))
       #f)
      ((equal? a c)
       (or (equal? b d)
           (ordered-ccw (reference-direction a) d b a)))
      ((equal? b d)
       (ordered-ccw (reference-direction b) c a b))
      ((equal? a d)
       (or (equal? b c)
           (ordered-ccw (reference-direction a) c b a)))
      ((equal? b c)
       (ordered-ccw (reference-direction b) d a b))
      (#t
       (error "vertex-crossing? called with 4 distinct vertices"))))
  result)

;; Return true if the segments S1 and S2 are crossing, even if they share a
;; common vector.
(: segment-or-vertex-crossing? (-> segment segment Boolean))
(define (segment-or-vertex-crossing? s1 s2)
  (define sign (segment-crossing-sign s1 s2))
  (if (zero? sign)
      (vertex-crossing? s1 s2)
      (> sign 0)))

;; The origin for ray casting when determining if a point is inside the loop.
;; The value is chosen to minimize the number of expensive triangle sign
;; operations, but otherwise it is arbitrary.  This is the same value as the
;; one used in the S2 library.
(: origin FlVector)
(define origin (flvector -0.0099994664350250197 0.0025924542609324121 0.99994664350250195))

;; Return a suitable origin location for ray casting tests for a loop defined
;; by SEGMENTS as well as a boolean value indicating if this origin is inside
;; or outside the loop.
(: get-suitable-origin (-> (Listof segment) (values Vector3 Boolean)))
(define (get-suitable-origin segments)
  (define-values (v0 v1 v2)
    (let ([s0 (car segments)]
          [s1 (cadr segments)])
      (values (segment-start s0) (segment-end s0) (segment-end s1))))
  (values
   origin
   (let ([v1-inside?
          (and (not (equal? v0 v1))
               (not (equal? v1 v2))
               (angle-contains-vertex? v0 v1 v2))]
         [contains-v1? (is-point-inside-loop? v1 segments origin #f)])
     (define is-inside? (not (equal? v1-inside? contains-v1?)))
     is-inside?)))

;; Return true if the point P is inside the loop formed by the list of
;; SEGMENTS.  ORIGIN is used as the "point at infinity" for the segment
;; crossing test and ORIGIN-INSIDE? specifies if this origin point is inside
;; or outside the loop (it will usually be outside).
(: is-point-inside-loop? (-> Vector3 (Listof segment) Vector3 Boolean Boolean))
(define (is-point-inside-loop? p segments origin is-origin-inside?)
  (let ([tseg : segment (make-segment origin p)])
    (let loop : Boolean
         ([crossing? : Boolean is-origin-inside?]
          [remaining : (Listof segment) segments])
      (if (null? remaining)
          crossing?
          (let ([x? (segment-or-vertex-crossing? tseg (car remaining))])
            (loop (assert (xor crossing? x?) boolean?)
                  (cdr remaining)))))))


;;........................................................... The Region ....

;; The types of cell - region containment relationships, as reported by
;; `relationship-to-cell`
(define-type Cell-Region-Relationship
  (U 'region-inside-cell                ; the region is inside the cell
     'region-contains-cell              ; the region contains the cell
     'region-intersects-cell            ; the region and cell intersect
     'none                              ; the region and cell are unrelated
     ))

;; The region represents a region on the earth surface, for which we'll do
;; geoid tiling.  Specific regions are being implemented below, but this is
;; the common interface that the tiling code uses.
(define-type Region%
  (Class
   ;; Return a list of geoids from where the tiling will start. Note that
   ;; there might be multiple starting geoids if a region as multiple
   ;; non-overlapping sub-regions
   [initial-geoids (-> Integer (Listof Integer))]
   ;; Determine the relationship between this region and a cell
   [relationship-to-cell (-> cell Cell-Region-Relationship)]))

(: region% Region%)
(define region%
  (class object%
    (super-new)
    (define/public (initial-geoids _level) null)
    (define/public (relationship-to-cell _cell) 'none)))


;;.................................................... The Spherical Cap ....

(define-type Sperical-Cap%
  (Class
   #:implements Region%
   (init-field
    [center Vector3]
    [chlen² Flonum])
   [get-center (-> Vector3)]
   [get-radius (-> Flonum)]
   [is-empty? (-> Boolean)]
   [is-full? (-> Boolean)]
   [intersects-cell?/no-vertices (-> cell Boolean)]))

(: sperical-cap% Sperical-Cap%)
(define sperical-cap%
  (class region%
    (init-field center chlen²)
    (super-new)

    (: the-complement (U #f (Instance Sperical-Cap%)))
    (define the-complement #f)

    (define/override (initial-geoids level)
      (if (is-empty?)
          '()
          (list (unit-vector->geoid (flvector-ref center 0)
                                    (flvector-ref center 1)
                                    (flvector-ref center 2)
                                    level))))

    (define/override (relationship-to-cell cell)
      (define inside-vertex-count
        (for/sum : Integer ([p : Vector3 (in-list (cell-corners cell))])
          (if (is-point-inside? p) 1 0)))
      (cond
        ((= inside-vertex-count 4)     ; all vertices are inside the cap
         (if (send (complement) intersects-cell?/no-vertices cell)
             'region-intersects-cell
             'region-contains-cell))
        ((> inside-vertex-count 0)      ; some vertices are inside
         'region-intersects-cell)
        (#t                             ; no vertices are inside the cap
         (cond ((intersects-cell?/no-vertices cell)
                'region-intersects-cell)
               ((point-inside-cell? center cell)
                'region-inside-cell)
               (#t
                'none)))))

    (define/public (get-center) center)
    (define/public (get-radius) (chlen²->radians chlen²))

    (define/public (is-empty?)
      ;; Note that a chord length of 0 is a cap containing only the center!
      (< chlen² 0.0))

    (define/public (is-full?)
      (>= chlen² 4.0))

    (: complement (-> (Instance Sperical-Cap%)))
    (define/private (complement)
      (cond ((is-empty?) the-full-cap)
            ((is-full?) the-empty-cap)
            (#t
             (unless the-complement
               (set! the-complement (new sperical-cap% [center (v3- center)] [chlen² (- 4.0 chlen²)])))
             ;; See https://racket.discourse.group/t/type-checking-class-instances/1044
             (define c the-complement)
             (if c c (error "sperical-cap%/complement: should not happen!")))))

    (: is-point-inside? (-> Vector3 Boolean))
    (define/private (is-point-inside? p)
      (or (is-full?)
          (<= (chlen²-from-points center p) chlen²)))

    (: sin2 (-> Flonum))
    (define/private (sin2)
      (* chlen² (- 1.0 (* 0.25 chlen²))))

    (define/public (intersects-cell?/no-vertices cell)
      (cond
        ((>= chlen² 2.0)
         ;; cap is larger than the hemisphere, and we already know the cell
         ;; vertices are not inside the cap
         #f)
        ((is-empty?)
         #f)
        ((point-inside-cell? center cell)
         #t)
        (#t
         (define sin2a (sin2))
         (let loop ([edges (cell-edges/raw cell)]
                    [corners (cell-corners cell)])
           (if (or (null? edges) (null? corners))
               #f
               (let* ([edge (car edges)]
                      [dot (dot-product center edge)])
                 (cond ((> dot 0)
                        (loop (cdr edges) (cdr corners)))
                       ((> (* dot dot) (* sin2a (v3len² edge)))
                        #f)
                       (#t
                        (let ([dir (cross-product edge center)])
                          (assert edges pair?)
                          (let ([c1 (car corners)]
                                [c2 (if (null? (cdr corners)) (car (cell-corners cell)) (cadr corners))])
                            (if (and (< (dot-product dir c1) 0.0) (> (dot-product dir c2) 0.0))
                                #t
                                (loop (cdr edges) (cdr corners)))))))))))))
    ))

(define the-empty-cap (new sperical-cap% [center (flvector 0.0 0.0 1.0)] [chlen² -1.0]))
(define the-full-cap (new sperical-cap% [center (flvector 0.0 0.0 1.0)] [chlen² 4.0]))

(: make-spherical-cap (-> Real Real Real (Instance Sperical-Cap%)))
(define (make-spherical-cap lat lng radius)
  (define-values (x y z) (lat-lng->unit-vector lat lng))
  (define a (/ (real->double-flonum radius) earth-radius))
  (new sperical-cap%
       [center (flvector x y z)]
       [chlen² (radians->chlen² a)]))


;;.................................................... The Open Polyline ....

(define-type Open-Polyline%
  (Class
   #:implements Region%
   (init-field
    [segments (Listof segment)])))

(: open-polyline% Open-Polyline%)
(define open-polyline%
  (class region%
    (init-field segments)
    (super-new)

   (define/override (initial-geoids level)
     (if (null? segments)
         null
         (let ([first-segment (car segments)])
           (list (segment->start-geoid first-segment level)))))

   (define/override (relationship-to-cell cell)
     ;; Note that a cell cannot be inside an open poliline, so
     ;; 'region-contains-cell will never be returned here.
     (cond ((ormap (lambda ([s : segment]) : Boolean
                     (segment-intersects-cell?/strict s cell))
                   segments)
            'region-intersects-cell)
           ;; NOTE: we assume connected segments, so we check points once
           ;; only.
           ((and (not (null? segments))
                 (point-inside-cell? (segment-start (car segments)) cell)
                 (andmap (lambda ([s : segment]) : Boolean
                           (point-inside-cell? (segment-end s) cell))
                         segments))
            'region-inside-cell)
           (#t
            'none)))

    ))


(: make-open-polyline (-> (Listof (Vector Real Real)) (Instance Open-Polyline%)))
(define (make-open-polyline track)
  (new open-polyline% [segments (track->segments track)]))


;;.................................................. The Closed Polyline ....

(define-type Closed-Polyline%
  (Class
   #:implements Region%
   (init-field
    [segments (Listof segment)])))

(: closed-polyline% Closed-Polyline%)
(define closed-polyline%
  (class region%
    (init-field segments)

    (: sentry-level Integer)
    (define sentry-level -1)

    (: origin Vector3)
    (: is-origin-inside? Boolean)
    (define-values (origin is-origin-inside?)
      (get-suitable-origin segments))

    (define/override (initial-geoids level)
      (if (null? segments)
          null
          (let ([first-segment (car segments)])
            (list (segment->start-geoid first-segment level)))))

    (define/override (relationship-to-cell cell)
      (cond ((loop-intersects-cell? cell)
             'region-intersects-cell)
            ((loop-contains-some-cell-corners? cell)
             (when (> (cell-level cell) sentry-level)
               (set! sentry-level (cell-level cell)))
             'region-contains-cell)
            ;; A cell can contain this region only if its level is larger than
            ;; the level of the largest cell we found to fit inside the
            ;; region.
            ((and (> (cell-level cell) sentry-level)
                  (loop-inside-cell? cell))
             'region-inside-cell)
            (#t
             'none)))

    (: any-segment-intersects-cell? (-> (Listof segment) cell Boolean))
    (define/private (any-segment-intersects-cell? segments cell)
      (ormap (lambda ([s : segment]) : Boolean
               (segment-intersects-cell?/strict s cell)) segments))

    (: loop-intersects-cell? (-> cell Boolean))
    (define/private (loop-intersects-cell? cell)
      (any-segment-intersects-cell? segments cell))

    (: loop-inside-cell? (-> cell Boolean))
    (define/private (loop-inside-cell? cell)
      (and (not (null? segments))
           (point-inside-cell? (segment-start (car segments)) cell)
           (andmap (lambda ([s : segment]) : Boolean
                     (point-inside-cell? (segment-end s) cell))
                   segments)))

    (: loop-contains-some-cell-corners? (-> cell Boolean))
    (define/private (loop-contains-some-cell-corners? cell)
      (ormap point-inside-loop? (cell-corners cell)))

    (: point-inside-loop? (-> Vector3 Boolean))
    (define (point-inside-loop? p)
      (is-point-inside-loop? p segments origin is-origin-inside?))

    (super-new)

    ))


;;.......................................... The Indexed Closed Polyline ....

(: closed-polyline-indexed% Closed-Polyline%)
(define closed-polyline-indexed%
  (class region%
    (init-field segments)

    (: sentry-level Integer)
    (define sentry-level -1)

    (: key-level Integer)
    (define key-level 20)

    (: segment-index (HashTable Integer (Listof segment)))
    (define segment-index (make-hash))
    (: index-cells (HashTable Integer cell))
    (define index-cells (make-hash))

    (: origin Vector3)
    (: is-origin-inside? Boolean)
    (define-values (origin is-origin-inside?)
      (get-suitable-origin segments))

    (: build-index (-> Void))
    (define/private (build-index)
      (hash-clear! segment-index)
      (hash-clear! index-cells)
      (for ([s : segment (in-list segments)])
        (define covering (geoid-covering-for-segment s key-level))
        (for ([g : Integer (in-list covering)])
          (hash-update!
           segment-index
           g
           (lambda ([existing : (Listof segment)]) : (Listof segment) (cons s existing))
           (lambda () null))
          (hash-set! index-cells g (make-cell g)))))

    (define/override (initial-geoids level)
      (if (null? segments)
          null
          (let ([first-segment (car segments)])
            (list (segment->start-geoid first-segment level)))))

    (define/override (relationship-to-cell cell)
      (cond ((loop-intersects-cell? cell)
             'region-intersects-cell)
            ((loop-contains-some-cell-corners? cell)
             (when (> (cell-level cell) sentry-level)
               (set! sentry-level (cell-level cell)))
             'region-contains-cell)
            ;; A cell can contain this region only if its level is larger than
            ;; the level of the largest cell we found to fit inside the
            ;; region.
            ((and (> (cell-level cell) sentry-level)
                  (loop-inside-cell? cell))
             'region-inside-cell)
            (#t
             'none)))

    (: any-segment-intersects-cell? (-> (Listof segment) cell Boolean))
    (define/private (any-segment-intersects-cell? segments cell)
      (ormap (lambda ([s : segment]) : Boolean
               (segment-intersects-cell?/strict s cell)) segments))

    (: loop-intersects-cell? (-> cell Boolean))
    (define/private (loop-intersects-cell? cell)
      (if (> (cell-level cell) key-level)
          (any-segment-intersects-cell? segments cell)
          (let ([key (enclosing-geoid (cell-id cell) key-level)])
            (define candidates (hash-ref segment-index key (lambda () null)))
            (if (= key-level (cell-level cell))
                (not (null? candidates)) ; easy win, to avoid repeat
                (any-segment-intersects-cell? candidates cell)))))

    (: loop-inside-cell? (-> cell Boolean))
    (define/private (loop-inside-cell? cell)
      (and (not (null? segments))
           (point-inside-cell? (segment-start (car segments)) cell)
           (andmap (lambda ([s : segment]) : Boolean
                     (point-inside-cell? (segment-end s) cell))
                   segments)))

    (: loop-contains-some-cell-corners? (-> cell Boolean))
    (define/private (loop-contains-some-cell-corners? cell)
      (ormap point-inside-loop? (cell-corners cell)))

    (: point-inside-loop? (-> Vector3 Boolean))
    (define (point-inside-loop? p)
      (let ([test-segment : segment (make-segment origin p)])
        (for/fold ([crossing? : Boolean is-origin-inside?]
                   [visited : (Setof segment) (seteq)]
                   #:result crossing?)
                  ([(geoid cell) (in-hash index-cells)])
          (if (segment-intersects-cell?/relaxed test-segment cell)
              (let loop : (values Boolean (Setof segment))
                   ([crossing? : Boolean crossing?]
                    [visited : (Setof segment) visited]
                    [remaining : (Listof segment) (hash-ref segment-index geoid (lambda () null))])
                (if (null? remaining)
                    (values crossing? visited)
                    (let ([candidate (car remaining)])
                      (if (set-member? visited candidate)
                          (loop crossing? visited (cdr remaining))
                          (let ([x? (segment-or-vertex-crossing? test-segment candidate)])
                            (loop (assert (xor crossing? x?) boolean?)
                                  (set-add visited candidate)
                                  (cdr remaining)))))))
              (values crossing? visited)))))

    (super-new)
    (build-index)

    ))

(: make-closed-polyline (->* ((Listof (Vector Real Real)))
                             (#:ccw? Boolean #:force-indexed? Boolean)
                             (Instance Closed-Polyline%)))
(define (make-closed-polyline track
                              #:ccw? [ccw? #t]
                              #:force-indexed? [force-indexed? #f])
  (when (< (length track) 3)
    (error "make-closed-polyline: track must contain at least 3 points"))
  (define segments
    (track->segments (if ccw? track (reverse track)) #t))
  ;; This can be very expensive to call...
  #;(unless (segments-valid-for-loop? segments)
    (error "make-closed-polyline: cannot handle self-intersecting track"))
  (if (or force-indexed? (> (length segments) 100))
      (new closed-polyline-indexed% [segments segments])
      (new closed-polyline% [segments segments])))


;;..................................................... The Union Region ....

(define-type Union-Region%
  (Class
   #:implements Region%
   (init-field
    [region-a (Instance Region%)]
    [region-b (Instance Region%)])))

(: union-region% Union-Region%)
(define union-region%
  (class object%
    (init-field region-a region-b)
    (super-new)

    (define/public (initial-geoids level)
      (remove-duplicates
       (append (send region-a initial-geoids level)
               (send region-b initial-geoids level))))

    (define/public (relationship-to-cell cell)
      (let ([ra (send region-a relationship-to-cell cell)])
        (case ra
          ((region-contains-cell)
           'region-contains-cell)
          ((region-intersects-cell)
           (let ([rb (send region-b relationship-to-cell cell)])
             (case rb
               ((region-contains-cell)
                'region-contains-cell)
               (else
                'region-intersects-cell))))
          ((region-inside-cell)
           (let ([rb (send region-b relationship-to-cell cell)])
             (case rb
               ((region-inside-cell)
                'region-inside-cell)
               (else
                'region-intersects-cell))))
          ((none)
           (let ([rb (send region-b relationship-to-cell cell)])
             ;; So region-a returned 'none on this cell and region-b didn't.
             ;; Swap them, so next time we ask region-b first.
             (unless (equal? rb 'none)
               (let ([tmp region-a])
                 (set! region-a region-b)
                 (set! region-b tmp)))

             (case rb
               ((region-inside-cell)
                'region-intersects-cell)
               (else
                rb)))))))
    ))

(: join-regions (->* () () #:rest (Instance Region%) (Instance Region%)))
(define (join-regions . regions)
  (if (null? regions)
      (new region%)                     ; this is an empty region
      (let loop : (Instance Region%)
           ([result : (Instance Region%) (car regions)]
            [remaining : (Listof (Instance Region%)) (cdr regions)])
        (if (null? remaining)
            result
            (loop (new union-region%
                       [region-a result]
                       [region-b (car remaining)])
                  (cdr remaining))))))


;;................................................. The Intersect Region ....

(define-type Intersect-Region%
  (Class
   #:implements Region%
   (init-field
    [region-a (Instance Region%)]
    [region-b (Instance Region%)])))

(: intersect-region% Intersect-Region%)
(define intersect-region%
  (class object%
    (init-field region-a region-b)
    (super-new)

    (define/public (initial-geoids level)
      ;; Intersection has to be reachable from Region a
      (send region-a initial-geoids level))

    (define/public (relationship-to-cell cell)
      (let ([ra (send region-a relationship-to-cell cell)])
        (if (equal? ra 'none)
            'none
            (let ([rb (send region-b relationship-to-cell cell)])
              (case ra
                ((region-contains-cell)
                 (case rb
                   ((region-contains-cell region-inside-cell)
                    'region-contains-cell)
                   ((region-intersects-cell)
                    'region-intersects-cell)
                   ((none)
                    ;; So region-a returned something on this cell and
                    ;; region-b didn't.  Swap them, so next time we ask
                    ;; region-b first, which would be faster
                    (let ([tmp region-a])
                      (set! region-a region-b)
                      (set! region-b tmp))
                    'none)))
                ((region-intersects-cell)
                 (case rb
                   ((region-contains-cell region-intersects-cell)
                    'region-intersects-cell)
                   ((region-inside-cell)
                    'region-inside-cell)
                   ((none)
                    ;; So region-a returned something on this cell and
                    ;; region-b didn't.  Swap them, so next time we ask
                    ;; region-b first, which would be faster
                    (let ([tmp region-a])
                      (set! region-a region-b)
                      (set! region-b tmp))
                    'none)))
                ((region-inside-cell)
                 (case rb
                   ((region-contains-cell region-intersects-cell)
                    'region-intersects-cell)
                   ((region-inside-cell)
                    'region-inside-cell)
                   ((none)
                    ;; So region-a returned something on this cell and
                    ;; region-b didn't.  Swap them, so next time we ask
                    ;; region-b first, which would be faster
                    (let ([tmp region-a])
                      (set! region-a region-b)
                      (set! region-b tmp))
                    'none)))
                ((none)
                 'none))))))
    ))

(: intersect-regions (->* () () #:rest (Instance Region%) (Instance Region%)))
(define (intersect-regions . regions)
  (if (null? regions)
      (new region%)                     ; this is an empty region
      (let loop : (Instance Region%)
           ([result : (Instance Region%) (car regions)]
            [remaining : (Listof (Instance Region%)) (cdr regions)])
        (if (null? remaining)
            result
            (loop (new intersect-region%
                       [region-a result]
                       [region-b (car remaining)])
                  (cdr remaining))))))


;;.................................................. The Subtract Region ....

(define-type Subtract-Region%
  (Class
   #:implements Region%
   (init-field
    [region-a (Instance Region%)]
    [region-b (Instance Region%)])))

(: subtract-region% Subtract-Region%)
(define subtract-region%
  (class object%
    (init-field region-a region-b)
    (super-new)

    (define/public (initial-geoids level)
      ;; Subtraction has to be reachable from Region a
      (send region-a initial-geoids level))

    (define/public (relationship-to-cell cell)
      (define rb (send region-b relationship-to-cell cell))
      (case rb
        ((region-contains-cell) 'none)
        ((none) (send region-a relationship-to-cell cell))
        ((region-intersects-cell region-inside-cell)
         (let ([ra (send region-a relationship-to-cell cell)])
           (case ra
             ((region-contains-cell) 'region-intersects-cell)
             (else ra))))))))

(: subtract-regions (-> (Instance Region%) (Instance Region%) (Instance Region%)))
(define (subtract-regions ra rb)
  (new subtract-region% [region-a ra] [region-b rb]))


;;............................................................... tiling ....

(: coarse-geoid-tiling-for-region
   (-> (Instance Region%) Integer
       (values (Listof Integer) (Listof Integer))))
(define (coarse-geoid-tiling-for-region region max-level)
  (let loop : (values (Listof Integer) (Listof Integer))
       ([contains : (Setof Integer) (set)]
        [intersects : (Setof Integer) (set)]
        [frontier : (Listof Integer) (send region initial-geoids max-level)]
        [visited : (Setof Integer) (set)])
    (if (null? frontier)
        (values (set->list contains) (set->list intersects))
        (let ([geoid (car frontier)]
              [frontier (cdr frontier)])
          (if (set-member? visited geoid)
              (loop contains intersects frontier visited)
              (case (send region relationship-to-cell (make-cell geoid))
                ((region-inside-cell)
                 ;; Found a geoid that contains the entire region, no point in
                 ;; looking for others...  Note that we put this in the
                 ;; intersects set, to be refined later...
                 (values null (list geoid)))
                ((region-contains-cell)
                 (loop (set-add contains geoid)
                       intersects
                       (append (adjacent-geoids geoid) frontier)
                       (set-add visited geoid)))
                ((region-intersects-cell)
                 (loop contains
                       (set-add intersects geoid)
                       (append (adjacent-geoids geoid) frontier)
                       (set-add visited geoid)))
                ((none)
                 (loop contains
                       intersects
                       frontier
                       (set-add visited geoid)))))))))

(: refine-geoid-tiling-for-region
   (-> (Instance Region%) Integer Integer
       (Listof Integer)))
(define (refine-geoid-tiling-for-region region geoid min-level)

  (: as-list (-> (U Integer (Listof Integer)) (Listof Integer)))
  (define (as-list item) (if (integer? item) (list item) item))

  (: refine (-> Integer Integer (U Integer (Listof Integer))))
  (define (refine geoid level)
    (if (<= level min-level)
        geoid                           ; reached min-level, no more splitting
        (let ([candidates (split-geoid geoid)])
          (: splits (Listof (U Integer (Listof Integer))))
          (define splits
            (map (lambda ([g : Integer]) : (U Integer (Listof Integer))
                   (case (send region relationship-to-cell (make-cell g))
                     ((region-inside-cell region-intersects-cell)
                      (refine g (sub1 level)))
                     ((region-contains-cell)
                      g)
                     ((none)
                      null)))
                 candidates))
          (if (andmap integer? splits)
              geoid ; no need of subsplitting, all subdivisions are inside the region
              (append* (map as-list splits))))))

  (as-list (refine geoid (geoid-level geoid))))

(: geoid-tiling-for-region (->* ((Instance Region%) Integer Integer)
                                (#:verbose? Boolean)
                                (Listof Integer)))
(define (geoid-tiling-for-region region min-level max-level #:verbose? [verbose? #f])
  (define start (current-inexact-milliseconds))
  (let-values ([(contains intersects)
                (coarse-geoid-tiling-for-region region max-level)])
    (define duration (- (current-inexact-milliseconds) start))
    (when verbose?
      (printf "*** Completed coarse tiling at level ~a in ~a seconds (contains = ~a, intersects = ~a)~%"
              max-level (~r (/ duration 1000.0) #:precision 2) (length contains) (length intersects)))
    (: refine (-> Integer (Listof Integer)))
    (define (refine g)
      (define start (current-inexact-milliseconds))
      (define result (refine-geoid-tiling-for-region region g min-level))
      (when verbose?
        (define duration (- (current-inexact-milliseconds) start))
        (printf "*** Refining #x~x took ~a seconds producing ~a geoids~%"
                g (~r (/ duration 1000.0) #:precision 2) (length result)))
      result)
    (append contains (append* (map refine intersects)))))


;;.................................................. guess-winding-order ....

;; Return the angle between the vectors (x0, y0) -- (x1, y1) and (x0, y0) --
;; (x2, y2)
(: subtended-angle (-> Flonum Flonum Flonum Flonum Flonum Flonum Flonum))
(define (subtended-angle x0 y0 x1 y1 x2 y2)
  (define s1x (- x1 x0))
  (define s1y (- y1 y0))
  (define s1len (assert (sqrt (+ (* s1x s1x) (* s1y s1y))) flonum?))
  (define s2x (- x2 x0))
  (define s2y (- y2 y0))
  (define s2len (assert (sqrt (+ (* s2x s2x) (* s2y s2y))) flonum?))
  (define dot-product (+ (* s1x s2x) (* s1y s2y)))
  (if (or (zero? dot-product) (zero? s1len) (zero? s1len))
      0.0
      (let ([angle (acos (min 1.0 (/ dot-product (* s1len s2len))))])
        (define cross-magnitude (- (* (- x1 x0) (- y2 y0))
                                   (* (- y1 y0) (- x2 x0))))
        (* angle (flsgn cross-magnitude)))))

;; Guess the winding order (clockwise or counter-clockwise) for the polygon
;; defined by TRACK, which is a list of latitude/longitude points. The guess
;; works for "usual" polygons which are much smaller than a hemisphere and
;; assumes that the "inside" is the smaller part.
;;
;; GeoJSON objects don't specify a winding order, while our closed polylines
;; expect the data to be in CCW order (and they can do the conversion
;; themselves), so this function can be used to determine the winding order of
;; polygons inside GeoJSON objects.
;;
;; NOTE: on a sphere, such as Earth, each loop of points defines two regions;
;; the simplest case to understand is a region defined by a loop around the
;; equator.  This can define either the northern or the southern hemisphere
;; and the difference between the two is defined by the winding order of the
;; points: the region which has the points in CCW order is the one defined by
;; the loop.
(: guess-winding-order (-> (Listof (Vector Real Real)) (U 'cw 'ccw)))
(define (guess-winding-order track)
  (define-values (lat0 lon0)
    (match-let ([(vector lat lon) (car track)])
      (values
       (real->double-flonum lat)
       (real->double-flonum lon))))
  (define winding-angle
    (let loop : Flonum
         ([winding-angle : Flonum 0.0]
          [x1 : Flonum lon0]
          [y1 : Flonum lat0]
          [remaining : (Listof (Vector Real Real)) (cdr track)])
      (if (null? remaining)
          (+ winding-angle (subtended-angle lon0 lat0 x1 y1 lon0 lat0))
          (match-let ([(vector y2^ x2^) (car remaining)])
            (define x2 (real->double-flonum x2^))
            (define y2 (real->double-flonum y2^))
            (loop (+ winding-angle (subtended-angle lon0 lat0 x1 y1 x2 y2))
                  x2
                  y2
                  (cdr remaining))))))
  (if (> winding-angle 0) 'ccw 'cw))


;;............................................. segments-valid-for-loop? ....

;; Return true if the list of segments SEGS form a valid loop, that is, they
;; don't intersect.  Having self intersecting segments means we cannot perform
;; tiling operations.
(: segments-valid-for-loop? (-> (Listof segment) Boolean))
(define (segments-valid-for-loop? segs)

  (: intersects? (-> segment (Listof segment) Boolean))
  (define (intersects? s0 segments)
    (ormap (lambda ([s1 : segment]) : Boolean
             (cond
               ((or (equal? (segment-end s0) (segment-start s1))
                    (equal? (segment-start s0) (segment-end s1)))
                ;; Segments are connected to each other
                #f)
               ((segment-or-vertex-crossing? s0 s1))
               (#t
                #f)))
           segments))

  (let loop ([remaining segs])
    (if (null? remaining)
        #t                          ; we are done, no self intersections found
        (if (intersects? (car remaining) (cdr remaining))
            #f
            (loop (cdr remaining))))))
