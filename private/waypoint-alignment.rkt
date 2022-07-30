;; waypoint-alignment.rkt -- determine how close are two GPS paths.
;;
;; This file is part of geoid -- work efficiently with geographic data
;; Copyright (c) 2021 Alex Harsányi <AlexHarsanyi@gmail.com>
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

#lang typed/racket/base
(require racket/flonum
         "geoid.rkt")

(provide (all-defined-out))

(: ->unit-vectors (-> (Vectorof Integer) FlVector))
(define (->unit-vectors geoids)
  (define u (make-flvector (* (vector-length geoids) 3)))
  (for ([g (in-vector geoids)]
        [index (in-naturals)])
    (define base (* index 3))
    (define-values (x y z) (geoid->unit-vector g))
    (flvector-set! u base x)
    (flvector-set! u (+ base 1) y)
    (flvector-set! u (+ base 2) z))
  u)

;; Return the cost of the optimal mapping between two paths of waypoints
;; (expressed as vectors of geoids).  The lower the cost, the more closely the
;; two paths match each other (even if they don't have the same number of
;; waypoints).  The cost should be zero for a path compared against itself,
;; but we have an error of about 0.04 (4cm) per waypoint.  The cost will also
;; be proportional to the number of waypoints in the path, so to determine if
;; the two paths are the same, the length of the path needs to be taken into
;; account.
;;
;; This algorithm runs in quadratic time, O(vector-length SOURCE) x
;; (vector-length TARGET), and also uses memory proportional to this quadratic
;; time.
;;
;; Straightforward dynamic time warping implementation, currently used to
;; compare it to `dtw/memory-efficient`, but could also be used to extract the
;; actual mapping between waypoints in the two paths.
;;
;; https://en.wikipedia.org/wiki/Dynamic_time_warping

(: dtw (-> (Vectorof Integer) (Vectorof Integer) Real))
(define (dtw source target)
  (define s (->unit-vectors source))
  (define t (->unit-vectors target))

  (define n (vector-length source))
  (define m (vector-length target))

  (define dtw (make-flvector (* m n) +inf.0))

  (: vref (-> Integer Integer Flonum))
  (define (vref i j) (flvector-ref dtw (+ (* i m) j)))
  (: vset! (-> Integer Integer Flonum Void))
  (define (vset! i j v)
    (flvector-set! dtw (+ (* i m) j) v))
  (: unit-distance (-> Integer Integer Flonum))
  (define (unit-distance i j)
    (define dot
      (+ (* (flvector-ref s i) (flvector-ref t j))
         (* (flvector-ref s (+ i 1)) (flvector-ref t (+ j 1)))
         (* (flvector-ref s (+ i 2)) (flvector-ref t (+ j 2)))))
    (acos (max -1.0 (min 1.0 dot))))

  (vset! 0 0 (unit-distance 0 0))

  (for ([i (in-range 1 n)])
    (for ([j (in-range 1 m)])
      (define cost (unit-distance (* i 3) (* j 3)))
      (define min-previous-cost (min (vref (sub1 i) j)
                                     (vref i (sub1 j))
                                     (vref (sub1 i) (sub1 j))))
      (vset! i j (+ cost min-previous-cost))))

  (* earth-radius (vref (sub1 n) (sub1 m))))

;; Same as `dtw`, but uses less memory (and not quadratic to the size of the
;; input), will also be slightly faster, but still runs in quadratic time.
(: dtw/memory-efficient (-> (Vectorof Integer) (Vectorof Integer) Real))
(define (dtw/memory-efficient source target)

  ;; Ensure that the shortest vector is in M (and T), so we have the shortest
  ;; scan lines (THIS-ROW and CURRENT-LINE)
  (define-values (n m s t)
    (let ([n (vector-length source)]
          [m (vector-length target)]
          [s (->unit-vectors source)]
          [t (->unit-vectors target)])
      (if (< n m)
          (values m n t s)
          (values n m s t))))

  ;; DTW only needs data from the previous row as it advances, so we don't
  ;; need a full MxN matrix, just two rows, the current one and the previous
  ;; one.  Here we create two rows and we'll swap between them as the
  ;; algorithm progresses, effectively using only two rows to fill in the full
  ;; MxN matrix.
  ;;
  ;; This strategy means we won't be able to extract the waypoint alignment
  ;; pairs themselves, just the cost, but often this is what is needed
  ;; (i.e. to determine how close two paths are to each other.
  (define this-row (make-flvector m +inf.0))
  (define previous-row (make-flvector m +inf.0))

  ;; Return the distance on the unit sphere between the geoids at positions I
  ;; and J (in the source and target arrays).  For speed we'll use the S and T
  ;; arrays of geoids converted to unit vectors.
  (: unit-distance (-> Integer Integer Flonum))
  (define (unit-distance i j)
    (define dot
      (+ (* (flvector-ref s i) (flvector-ref t j))
         (* (flvector-ref s (+ i 1)) (flvector-ref t (+ j 1)))
         (* (flvector-ref s (+ i 2)) (flvector-ref t (+ j 2)))))
    (acos (max -1.0 (min 1.0 dot))))

  ;; Fill the costs in the CURRENT row (at ROW-INDEX), using data from the
  ;; PREVIOUS row.  CURRENT and PREVIOUS are one of the THIS-ROW and
  ;; PREVIOUS-ROW vectors (and they swap around each iteration).
  (: fill-row-with-costs (-> Integer FlVector FlVector Void))
  (define (fill-row-with-costs row-index current previous)
    (for ([j (in-range 1 m)])
      (define cost (unit-distance (* row-index 3) (* j 3)))
      (define min-previous-cost (min (flvector-ref previous j)
                                     (flvector-ref current (sub1 j))
                                     (flvector-ref previous (sub1 j))))
      (flvector-set! current j (+ cost min-previous-cost))))

  ;; Bootstrap the algorithm -- setup the first element in row 0 to the
  ;; distance between the first geoids in each path and fill in row 1.
  (flvector-set! previous-row 0 (unit-distance 0 0))
  (fill-row-with-costs 1 this-row previous-row)
  (flvector-set! previous-row 0 +inf.0)

  ;; Fill in the remaining rows always swapping between the two rows and
  ;; returning the last one that we filled in.
  (define last-line
    (ann (let loop ([i : Integer 2]
                    [this-row : FlVector previous-row]
                    [previous-row : FlVector this-row])
           (if (< i n)
               (begin
                 (fill-row-with-costs i this-row previous-row)
                 (loop (add1 i) previous-row this-row))
               ;; Note that we swapped lines already, so we need to return the
               ;; previous one.
               previous-row))
         FlVector))
  (* earth-radius (flvector-ref last-line (sub1 m))))

(define waypoint-alignment-cost dtw/memory-efficient)
