#lang typed/racket

;; This file is part of geoid -- work efficiently with geographic data
;; Copyright (c) 2020 Alex Harsányi <AlexHarsanyi@gmail.com>
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

(provide (all-defined-out))


;;.................................................... The Hilbert Curve ....

;; For an explanation of the Hilbert curve, see
;; https://en.wikipedia.org/wiki/Hilbert_curve

;; Masks are generated using:
;; (for ([i (in-range 0 32)]) (printf "#x~x~%" (sub1 (expt 2 i))))

(: masks (Immutable-Vectorof Positive-Integer))
(define masks #(#x1 #x3 #x7 #xf #x1f #x3f #x7f #xff #x1ff #x3ff #x7ff
                #xfff #x1fff #x3fff #x7fff #xffff #x1ffff #x3ffff #x7ffff
                #xfffff #x1fffff #x3fffff #x7fffff #xffffff #x1ffffff
                #x3ffffff #x7ffffff #xfffffff #x1fffffff #x3fffffff
                #x7fffffff))

(: rotate (-> Integer
              Integer Integer
              Boolean Boolean
              (Values Integer Integer)))
(define (rotate level x y rx ry)
  (if ry
      (values x y)
      (if rx
          (let* ([mask (vector-ref masks level)]
                 [nx (- mask x)]
                 [ny (- mask y)])
            (when (or (< nx 0) (< ny 0))
              (error (format "Unexpected negative result: nx = ~a, ny = ~a, level = ~a, mask = ~a, x = ~a, y = ~a"
                             nx ny level mask x y)))
            (values ny nx))
          (values y x))))

(: xy->hilbert-distance (-> Integer
                            Integer Integer
                            Integer))
(define (xy->hilbert-distance level x y)
  (let loop ([l : Integer level]
             [x : Integer x]
             [y : Integer y])
    (if (< l 0)
        0
        (let ([rx (bitwise-bit-set? x l)]
              [ry (bitwise-bit-set? y l)])
          (define val (cond ((and rx ry) 2) (rx 3) (ry 1) (else 0)))
          (define-values (nx ny) (rotate level x y rx ry))
          (+ (arithmetic-shift val (* 2 l))
             (loop (sub1 l) nx ny))))))

(: hilbert-distance->xy (-> Integer
                            Integer
                            (Values Integer Integer)))
(define (hilbert-distance->xy level d)
  (let loop ([current : Integer 0]
             [x : Integer 0]
             [y : Integer 0])
    (if (> current level)
        (values x y)
        (let-values
            ([(rx ry)
              (case (bitwise-and (arithmetic-shift d (- (* current 2))) 3)
                ((3) (values #t #f))
                ((2) (values #t #t))
                ((1) (values #f #t))
                ((0) (values #f #f))
                (else (error "Unexpected (2)")))])
          (define bit (arithmetic-shift 1 current))
          (define ux (+ (if rx bit 0) x))
          (define uy (+ (if ry bit 0) y))
          (define-values (nx ny) (rotate current ux uy rx ry))
          (loop (add1 current) nx ny)))))


;;.................................................. The Geometric Plane ....

;; A geometric plane defined by 4 real numbers, see
;; https://en.wikipedia.org/wiki/Plane_(geometry)

(struct plane ([a : Real] [b : Real] [c : Real] [d : Real]) #:transparent)

;; Construct a plane from an origin point (ox, oy, oz) and a normal vector
;; (nx, ny, nz).  The normal vector is assumed (but not checked) to be of unit
;; length.
(: make-plane (-> Real Real Real Real Real Real plane))
(define (make-plane ox oy oz nx ny nz)
  (plane nx ny nz (- (+ (* ox nx) (* oy ny) (* oz nz)))))

;; Return the distance from origin to the plane P along the vector x y z --
;; the vector is assumed to be of unit length.  Returns +inf.0 if the vector
;; is parallel to the plane.  Can return a negative value if the plane is in
;; the opposite direction from where the vector is pointing.
(: distance-to-plane (-> Real Real Real plane Real))
(define (distance-to-plane x y z p)
  (match-define (plane a b c d) p)
  (define div (+ (* x a) (* y b) (* z c)))
  (if (zero? div) +inf.0 (/ (- d) div)))


;;............................................................. The Face ....

;; A face represents one of our projection planes.  It has an identifier, plus
;; functions to compute the distance to the plane and to project / un-project
;; values to the plane.
(struct face ([id : Integer]
              [xyzd->uv : (-> Real Real Real Real (Values Real Real))]
              [uv->xyz : (-> Real Real (Values Real Real Real))]
              [distance : (-> Real Real Real Real)])
  #:transparent)

(: make-face (-> Integer
                 (-> Real Real Real Real (Values Real Real))
                 (-> Real Real (Values Real Real Real))
                 (-> Real Real Real Real)
                 face))
(define (make-face id xyzd->uv uv->xyz distance)
  (face id xyzd->uv uv->xyz distance))

;; NOTE: plane faces are numbered such that they are adjacent to each other
;; and the hilbert points continue from one face to the next: Front -> Right
;; -> Top -> Back -> Left -> Bottom, as shown on this page:
;;
;; http://s2geometry.io/devguide/s2cell_hierarchy

(define front-face
  (make-face
   0
   (lambda (x y z d)
     (values (* 1/2 (+ 1 (* y d))) (* 1/2 (+ 1 (* z d)))))
   (lambda (u v)
     (define y (- (* 2 u) 1))
     (define z (- (* 2 v) 1))
     (define len (cast (sqrt (+ 1 (* y y) (* z z))) Real))
     (values (/ 1 len) (/ y len) (/ z len)))
   (let ([p (make-plane 1 0 0 1 0 0)])
     (lambda (x y z) (distance-to-plane x y z p)))))

(define right-face
  (make-face
   1
   (lambda (x y z d)
     (values (* 1/2 (+ 1 (* z d))) (* 1/2 (- 1 (* x d)))))
   (lambda (u v)
     (define z (- (* 2 u) 1))
     (define x (- 1 (* 2 v)))
     (define len (cast (sqrt (+ 1 (* x x) (* z z))) Real))
     (values (/ x len) (/ 1 len) (/ z len)))
   (let ([p (make-plane 0 1 0 0 1 0)])
     (lambda (x y z) (distance-to-plane x y z p)))))

(define top-face
  (make-face
   2
   (lambda (x y z d)
     (values (* 1/2 (- 1 (* x d))) (* 1/2 (- 1 (* y d)))))
   (lambda (u v)
     (define x (- 1 (* 2 u)))
     (define y (- 1 (* 2 v)))
     (define len (cast (sqrt (+ (* x x) (* y y) 1)) Real))
     (values (/ x len) (/ y len) (/ 1 len)))
   (let ([p (make-plane 0 0 1 0 0 1)])
     (lambda (x y z) (distance-to-plane x y z p)))))

(define back-face
  (make-face
   3
   (lambda (x y z d)
     (values (* 1/2 (- 1 (* y d))) (* 1/2 (- 1 (* z d)))))
   (lambda (u v)
     (define y (- 1 (* 2 u)))
     (define z (- 1 (* 2 v)))
     (define len (cast (sqrt (+ 1 (* y y) (* z z))) Real))
     (values (/ -1 len) (/ y len) (/ z len)))
   (let ([p (make-plane -1 0 0 -1 0 0)])
     (lambda (x y z) (distance-to-plane x y z p)))))

(define left-face
  (make-face
   4
   (lambda (x y z d)
     (values (* 1/2 (- 1 (* z d))) (* 1/2 (+ 1 (* x d)))))
   (lambda (u v)
     (define z (- 1 (* 2 u)))
     (define x (- (* 2 v) 1))
     (define len (cast (sqrt (+ (* x x) 1 (* z z))) Real))
     (values (/ x len) (/ -1 len) (/ z len)))
   (let ([p (make-plane 0 -1 0 0 -1 0)])
     (lambda (x y z) (distance-to-plane x y z p)))))

(define bottom-face
  (make-face
   5
   (lambda (x y z d)
     (values (* 1/2 (+ 1 (* x d))) (* 1/2 (+ 1 (* y d)))))
   (lambda (u v)
     (define x (- (* 2 u) 1))
     (define y (- (* 2 v) 1))
     (define len (cast (sqrt (+ (* x x) (* y y) 1)) Real))
     (values (/ x len) (/ y len) (/ -1 len)))
   (let ([p (make-plane 0 0 -1 0 0 -1)])
     (lambda (x y z) (distance-to-plane x y z p)))))

;; NOTE: must be in face id order
(: all-faces (Listof face))
(define all-faces
  (list front-face right-face top-face back-face left-face bottom-face))

(: find-face (-> Real Real Real (values Real face)))
(define (find-face x y z)
  (define-values (d p)
    (for/fold ([d : Real +inf.0]
               [c : (U #f face) #f])
              ([p (in-list all-faces)])
      (define e ((face-distance p) x y z))
      (if (and (> e 0) (< e d))
          (values e p)
          (values d c))))
  (values d (cast p face)))


;;........................................................... level-info ....

;; Maximum encoding level.  0 is the level with the highest detail (a leaf
;; geoid)
(define max-level 30)                   ; 0 - 29
(define full-mask #xffffffffffffffff)

;; Store some information about each geoid level, which might be expensive to
;; compute.
(struct level-info
  ([id : Integer] ; the level number
   [mask : Integer] ; the mask value to clear the bits of the sentinel part
   [sentinel : Integer] ; the sentinel mask, identifying the level
   [max-coord : Integer] ; maximum x,y encoding coordinates at this level
   [epsilon : Real]) ; half square width at this level, relative to the unit face
  #:transparent)

(define level-information
  (for/vector : (Vectorof level-info)
      #:length (add1 max-level)
      ([id : Integer (in-range (add1 max-level))])
      (define flag-bit (* 2 id))
      (define sentinel (arithmetic-shift 1 flag-bit))
      (define mask
        (cast (sub1 (arithmetic-shift sentinel 1)) Integer))
      (define max-coord
        (arithmetic-shift 1 (- max-level id)))
      (define epsilon
        (if (> max-coord 1)
            (exact->inexact (* 0.5 (/ 1 (sub1 max-coord))))
            0))
      (level-info id mask sentinel max-coord epsilon)))


;;...................................................... The Unit Vector ....

;; Convert a latitude / longitude pair to a vector of unit length on the
;; unit-sphere.
(: lat-lng->unit-vector (-> Real Real (Values Real Real Real)))
(define (lat-lng->unit-vector lat lon)
  (define Θ (degrees->radians (- 90.0 lat)))
  (define Φ (degrees->radians lon))
  (define sin-Θ (sin Θ))
  (define cos-Θ (cos Θ))
  (define sin-Φ (sin Φ))
  (define cos-Φ (cos Φ))
  (values (* sin-Θ cos-Φ) (* sin-Θ sin-Φ) cos-Θ))

;; Convert a unit vector on the unit sphere back into the latitude / longitude
;; pair.
(: unit-vector->lat-lng (-> Real Real Real (Values Real Real)))
(define (unit-vector->lat-lng x y z)
  (define r (sqrt (+ (* x x) (* y y) (* z z))))
  (define cos-Θ (/ z r))
  (define Θ (cast (acos cos-Θ) Real))
  (define sin-Θ (sin Θ))
  (define cos-Φ (/ x sin-Θ))
  (define sin-Φ (/ y sin-Θ))
  (: Φ Number)               ; prevent TR from optimizing the real? call below
  (define Φ (acos cos-Φ))
  ;;(define Φ (atan (/ y x)))
  (define longitude (radians->degrees (if (real? Φ) Φ (real-part Φ))))
  ;; Note that we need bot the sine and cosine of Φ, to determine the correct
  ;; quadrant.
  (values (- 90.0 (radians->degrees Θ)) (if (< sin-Φ 0) (- longitude) longitude)))

;; Encode a latitude/longitude pair info the components of a geoid: the face,
;; the x, y coordinates on the plane for Hilbert distance computation and the
;; level.
(: encode (-> Real Real Integer
              (Values Integer
                      Integer
                      Integer
                      Integer)))
(define (encode lat lon level)
  (define-values (x y z) (lat-lng->unit-vector lat lon))

  (define-values (d f) (find-face x y z))
  (define-values (u v) ((face-xyzd->uv f) x y z d))
  ;; (printf "x = ~a, y = ~a, z = ~a; u = ~a, v = ~a~%" x y z u v)

  (match-define (level-info id mask sentinel max-coord epsilon)
    (vector-ref level-information level))
  (define semicircle-mask (sub1 max-coord))

  (define iu (max 0 (min semicircle-mask (exact-truncate (* u semicircle-mask)))))
  (define iv (max 0 (min semicircle-mask (exact-truncate (* v semicircle-mask)))))
  (values (face-id f)
          (cast iu Integer)
          (cast iv Integer)
          level))

;; Decode a latidude/longitude pair from the components of a geoid: the face,
;; the x, y integer coordinates and the level.
(: decode (-> Integer Integer Integer Integer
              (Values Real Real)))
(define (decode face ix iy level)
  (when (> face 5)
    ; there are only 6 faces, from 0 to 5, but the 3 bits used to encode them
    ; can also represent 6 and 7.  Note that (pack 6 0 0) is used as the
    ; sentinel geoh id (see below)
    (error (format "Bad face id: ~a" face)))

  (match-define (level-info id ask sentinel max-coord epsilon)
    (vector-ref level-information level))
  (define semicircle-mask (sub1 max-coord))

  (: u Real)
  (define u (+ epsilon (exact->inexact (/ ix semicircle-mask))))
  (: v Real)
  (define v (+ epsilon (exact->inexact (/ iy semicircle-mask))))
  (define f (list-ref all-faces face))
  (define-values (x y z) ((face-uv->xyz f) u v))
  ;; (printf "x = ~a, y = ~a, z = ~a; u = ~a, v = ~a~%" x y z u v)
  (unit-vector->lat-lng x y z))

;; Pack the geoid components into a single 64 bit integer.
(: pack (-> Integer Integer Integer Integer
            Integer))
(define (pack face ix iy level)
  (define sentinel-position (* 2 level))
  (bitwise-ior
   (arithmetic-shift face (add1 (* 2 max-level)))
   (arithmetic-shift
    (xy->hilbert-distance (sub1 (- max-level level)) ix iy)
    (add1 sentinel-position))
   (arithmetic-shift 1 sentinel-position)))

;; Unpack a 64bit integer representing a geoid into its components: the face,
;; the x, y integer coordinates on the plane and the level.
(: unpack (-> Integer
              (Values Integer Integer Integer Integer)))
(define (unpack d)
  (define face
    (bitwise-and (arithmetic-shift d (- (add1 (* 2 max-level)))) #x7))

  (define linfo
    (for/or : (U #f level-info)
        ([linfo (in-vector level-information)]
         #:when (= (bitwise-and d (level-info-mask linfo))
                   (level-info-sentinel linfo)))
      linfo))

  (if linfo
      (let ([hd (arithmetic-shift d (- (add1 (* 2 (level-info-id linfo)))))])
        ;; NOTE that hilbert-distance->xy is not confused by the high bits
        ;; that are set in d for the plane face , so we don't need to clear
        ;; them
        (define-values (ix iy)
          (hilbert-distance->xy
           (cast (sub1 (- max-level (level-info-id linfo))) Integer)
           hd))
        (values face ix iy (level-info-id linfo)))
      (error (format "unpack: cannot find level information for ~a" d))))


;;.................................................................. API ....

(: the-first-geoid Integer)
(define the-first-geoid (pack 0 0 0 0))
(: the-last-geoid Integer)
(define the-last-geoid
  (let ([level-zero (vector-ref level-information 0)])
    (pack 5 (sub1 (level-info-max-coord level-zero)) 0 0)))
(: the-sentinel-geoid Integer)
(define the-sentinel-geoid (add1 the-last-geoid))

(: first-valid-geoid (-> Integer))
(define (first-valid-geoid) the-first-geoid)

(: last-valid-geoid (-> Integer))
(define (last-valid-geoid) the-last-geoid)

(: sentinel-geoid (-> Integer))
(define (sentinel-geoid) the-sentinel-geoid)

(: valid-geoid? (-> Integer Boolean))
(define (valid-geoid? geoid)
  (and (>= geoid the-first-geoid) (<= geoid the-last-geoid)))

;; Produce a random GEOID at the specified LEVEL.  This is used by the test
;; suite.
(: random-geoid (-> Integer Integer))
(define (random-geoid level)
  (define linfo (vector-ref level-information level))
  (define face (random 6))
  (define ix (random (level-info-max-coord linfo)))
  (define iy (random (level-info-max-coord linfo)))
  (pack face ix iy level))

(: geoid-level (-> Integer Integer))
(define (geoid-level geoid)
  (unless (valid-geoid? geoid)
    (error (format "geoid-level: invalid geoid: ~a" geoid)))
  (define linfo
    (for/or : (U #f level-info)
        ([linfo (in-vector level-information)]
         #:when (= (bitwise-and geoid (level-info-mask linfo))
                   (level-info-sentinel linfo)))
      linfo))
  (if linfo
      (level-info-id linfo)
      (error (format "geoid-level: cannot find level information for ~a" geoid))))

;; Equivalent to (= (geoid-level geoid) 1), but faster
(: leaf-geoid? (-> Integer Boolean))
(define (leaf-geoid? geoid)
  ;; Leaf GEOIDS have 1 as their least significant bit, so they are always
  ;; oddd.
  (and (valid-geoid? geoid) (odd? geoid)))

(: geoid-face (-> Integer Integer))
(define (geoid-face geoid)
  (bitwise-and (arithmetic-shift geoid (- (add1 (* 2 max-level)))) #x7))

;; Return the integer value that you can add to a geoid to obtain the next
;; geoid at the same level.
(: geoid-stride (-> Integer Integer))
(define (geoid-stride geoid)
  (define linfo (vector-ref level-information (geoid-level geoid)))
  (add1 (level-info-mask linfo)))

(: enclosing-geoid (-> Integer Integer Integer))
(define (enclosing-geoid geoid new-level)
  (unless (valid-geoid? geoid)
    (raise-argument-error 'geoid "valid-geoid?" geoid))
  (when (or (< new-level 0) (>= new-level max-level))
    (raise-argument-error 'new-level (format "(between/c 0 ~a)" max-level) new-level))
  (let ([level (geoid-level geoid)])
    (cond ((= level new-level) geoid)
          ((>= level new-level)
           (raise-arguments-error
            'enclosing-geoid
            "input geoid at higher level"
            "geoid" geoid
            "geoid-level" level
            "new-level" new-level))
          (#t
           (let ([linfo (vector-ref level-information new-level)])
             (bitwise-ior (bitwise-and (bitwise-xor full-mask (level-info-mask linfo))
                                       geoid)
                          (level-info-sentinel linfo)))))))

(: lat-lng->geoid (-> Real Real Integer))
(define (lat-lng->geoid lat lon)
  (define-values (face ix iy level) (encode lat lon 0))
  (pack face ix iy level))

(: geoid->lat-lng (-> Integer (Values Real Real)))
(define (geoid->lat-lng id)
  (unless (valid-geoid? id)
    (raise-argument-error 'geoid "valid-geoid?" id))
  (define-values (face ix iy level) (unpack id))
  (decode face ix iy level))

;; Return the four leaf (level 0) geoids representing the corner of the geoid
;; ID.  Note that, while the GEOIDs are in CCW order, they are not necessarily
;; start from top-left and they are not in the order of their integer value.
(: leaf-corners (-> Integer (List Integer Integer Integer Integer)))
(define (leaf-corners id)
  (unless (valid-geoid? id)
    (raise-argument-error 'geoid "valid-geoid?" id))
  (define-values (face ix iy level) (unpack id))
  (define min-ix (arithmetic-shift ix level))
  (define min-iy (arithmetic-shift iy level))
  (define max-coord (sub1 (arithmetic-shift 1 level)))
  (list
   (pack face min-ix min-iy 0)
   (pack face min-ix (+ min-iy max-coord) 0)
   (pack face (+ min-ix max-coord) (+ min-iy max-coord) 0)
   (pack face (+ min-ix max-coord) min-iy 0)))

;; Split a geoid into the four geoids one level down.
(: split-geoid (-> Integer (List Integer Integer Integer Integer)))
(define (split-geoid id)
  (unless (valid-geoid? id)
    (raise-argument-error 'geoid "valid-geoid?" id))
  (when (leaf-geoid? id)
    (raise-argument-error 'geoid "(not leaf-geoid?)" id))
  (define-values (face ix iy level) (unpack id))
  (define min-ix (arithmetic-shift ix 1))
  (define min-iy (arithmetic-shift iy 1))
  (define new-level (sub1 level))
  (list
   (pack face min-ix min-iy new-level)
   (pack face min-ix (+ min-iy 1) new-level)
   (pack face (+ min-ix 1) (+ min-iy 1) new-level)
   (pack face (+ min-ix 1) min-iy new-level)))

;; Create an outline of the geoid ID as a list of leaf geoids going around the
;; edge in CCW direction.  As with `leaf-corners`, the first leaf geoid might
;; not necessarily be in the top left corner.  #:steps defines the number of
;; intermediate geoids to put around each side, wile #:closed?  defines
;; whether to add the first geoid as the last element in the loop.
(: leaf-outline (->* (Integer) (#:steps Integer #:closed? Boolean)
                     (Listof Integer)))
(define (leaf-outline id #:steps (steps 10) #:closed? (closed? #t))
  (unless (valid-geoid? id)
    (raise-argument-error 'geoid "valid-geoid?" id))
  (define-values (face ix iy level) (unpack id))
  (define min-ix (arithmetic-shift ix level))
  (define min-iy (arithmetic-shift iy level))
  (define max-coord (sub1 (arithmetic-shift 1 level)))
  (define stride (max 1 (exact-round (/ max-coord steps))))
  (append
   (for/list : (Listof Integer) ([delta (in-range 0 max-coord stride)])
     (pack face min-ix (+ min-iy delta) 0))
   (for/list : (Listof Integer) ([delta (in-range 0 max-coord stride)])
     (pack face (+ min-ix delta) (+ min-iy max-coord) 0))
   (for/list : (Listof Integer) ([delta (in-range max-coord 0 (- stride))])
     (pack face (+ min-ix max-coord) (+ min-iy delta) 0))
   (for/list : (Listof Integer) ([delta (in-range max-coord 0 (- stride))])
     (pack face (+ min-ix delta) min-iy 0))
   (if closed? (list (pack face min-ix min-iy 0)) '())))

;; Return the start and end leaf geoids betwen which all geoids inside this
;; one are located.
(: leaf-span (-> Integer (Values Integer Integer)))
(define (leaf-span id)
  (unless (valid-geoid? id)
    (raise-argument-error 'geoid "valid-geoid?" id))
  (define level (geoid-level id))
  (define base (- id (arithmetic-shift 1 (* 2 level)))) ; clear the level mask
  (values (+ base 1)                                    ; add the level 0 mask
          ;; Max value is either the sentinel-geoid, or the value which has
          ;; all lower bits set to 1 (including the last one which is the
          ;; level 0 mask)
          (min the-sentinel-geoid
               (+ 2 base (sub1 (arithmetic-shift 1 (add1 (* 2 level))))))))

(: contains-geoid? (-> Integer Integer Boolean))
(define (contains-geoid? this-geoid other-geoid)
  (unless (valid-geoid? this-geoid)
    (raise-argument-error 'this-geoid "valid-geoid?" this-geoid))
  (unless (valid-geoid? other-geoid)
    (raise-argument-error 'other-geoid "valid-geoid?" other-geoid))
  (define-values (start end) (leaf-span this-geoid))
  (and (>= other-geoid start) (< other-geoid end)))

;; Return the latitude/longitude rectangle for a GEOID.  Note that a geoid
;; does not correspond 1:1 to such a rectangle, and the returned rectangle
;; will contain more points than just this geoid.
(: lat-lng-rect (-> Integer (Values Real Real Real Real)))
(define (lat-lng-rect geoid)
  (unless (valid-geoid? geoid)
    (raise-argument-error 'geoid "valid-geoid?" geoid))

  ;; TODO: implement this later: the faces attain their min/max latitude and
  ;; longitude at their middle width and height.  Also top and bottom faces
  ;; include the poles and need special treatment.
  (when (= max-level (geoid-level geoid))
    (error "Cannot determine the lat-lng-rect for a face"))

  ;; NOTE: for non-face geoids, the max/min latitude and longitude is always
  ;; attained at the corners, but which ones depends on the face and the first
  ;; quadrant of the face -- this can be easily visualized if the few toplevel
  ;; geoids are drawn on a map using `leaf-outline`.
  ;;
  ;; We could save a few calculations by determining the face and quadrant and
  ;; only calculate latitude and longitude for the relevant corners
  ;; separately, but we keep things simple and just get the corners and find
  ;; the min/max from among the values.
  (match-define (list c0 c1 c2 c3) (leaf-corners geoid))
  (define-values (lat0 lng0) (geoid->lat-lng c0))
  (define-values (lat1 lng1) (geoid->lat-lng c1))
  (define-values (lat2 lng2) (geoid->lat-lng c2))
  (define-values (lat3 lng3) (geoid->lat-lng c3))
  ;; Enlarge the bounding box, by a small amount, to ensure the corners are
  ;; inside.
  (define e (* 2 (level-info-epsilon (vector-ref level-information 0))))
  (values
   (- (min lat0 lat1 lat2 lat3) e)
   (- (min lng0 lng1 lng2 lng3) e)
   (+ (max lat0 lat1 lat2 lat3) e)
   (+ (max lng0 lng1 lng2 lng3) e)))
