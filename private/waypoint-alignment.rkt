;; waypoint-alignment.rkt -- determine how close are two GPS paths.
;;
;; This file is part of geoid -- work efficiently with geographic data
;; Copyright (c) 2021, 2025 Alex Harsányi <AlexHarsanyi@gmail.com>
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
         racket/math
         "geoid.rkt")

(provide dtw
         dtw/memory-efficient
         dtw/window
         waypoint-alignment-cost)

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

;; Prepare the source and target vectors for running the DTW algorithms.
;; Given two geoid vectors, path1 path2, returns 5 values: the length of the
;; source, length of the target, and a unit-distance function which calculates
;; the angular distance between to unit vectors in the source and target
;; vectors.
;;
;; NOTE: source and target are always chosen such that source is the longest
;; vector between path1 and path2
(: prepare-source-and-target
   (-> (Vectorof Integer) (Vectorof Integer)
       (values Integer
               Integer
               (-> Integer Integer Flonum))))
(define (prepare-source-and-target path1 path2)

  ;; Ensure that the shortest vector is the target, in M (and T), so we have
  ;; the shortest scan lines (THIS-ROW and CURRENT-LINE)
  (define-values (n m s t)
    (let ([n (vector-length path1)]
          [m (vector-length path2)]
          [s (->unit-vectors path1)]
          [t (->unit-vectors path2)])
      (if (< n m)
          (values m n t s)
          (values n m s t))))

  ;; Return the distance on the unit sphere between the geoids at positions I
  ;; and J (in the source and target arrays).  For speed we'll use the S and T
  ;; arrays of geoids converted to unit vectors.
  (: unit-distance (-> Integer Integer Flonum))
  (define (unit-distance row col)
    (define i (* row 3))
    (define j (* col 3))
    (define dot
      (+ (* (flvector-ref s i) (flvector-ref t j))
         (* (flvector-ref s (+ i 1)) (flvector-ref t (+ j 1)))
         (* (flvector-ref s (+ i 2)) (flvector-ref t (+ j 2)))))
    (acos (max -1.0 (min 1.0 dot))))

  (values n m unit-distance))


;; Calculate window size for dtw/window based on the source and target vector
;; lengths.
(: get-window-size (-> Integer Integer Integer))
(define (get-window-size source-size target-size)
  ;; Input data should have been prepared by `prepare-source-and-target'
  (assert (and (>= source-size target-size) (>= target-size 3)))
  (define min-window-size 50)
  ;; A window size 5% of the target length seems to work for the data I have.
  ;; If the window is too small, we might obtain a higher cost than the
  ;; regular dtw cost, sometimes much higher.
  (define window-size (exact-round (* 0.05 target-size)))
  (min target-size (max min-window-size window-size)))

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
  (define-values (n m unit-distance)
    (prepare-source-and-target source target))

  (define dtw (make-flvector (* m n) +inf.0))

  (: vref (-> Integer Integer Flonum))
  (define (vref i j) (flvector-ref dtw (+ (* i m) j)))
  (: vset! (-> Integer Integer Flonum Void))
  (define (vset! i j v)
    (flvector-set! dtw (+ (* i m) j) v))

  ;; First element has nothing to look back to...
  (vset! 0 0 (unit-distance 0 0))
  ;; Rest of the elements in the first row also look back on the previous
  ;; element...
  (for ([j (in-range 1 m)])
    (define cost (unit-distance 0 j))
    (define min-previous-cost (vref 0 (sub1 j)))
    (vset! 0 j (+ cost min-previous-cost)))

  (for ([row-index (in-range 1 n)])

    ;; First element in each row only looks at the item above...
    (define cost (unit-distance row-index 0))
    (define min-previous-cost (vref (sub1 row-index) 0))
    (vset! row-index 0 (+ cost min-previous-cost))

    (for ([j (in-range 1 m)])
      (define cost (unit-distance row-index j))
      (define min-previous-cost
        (min (vref (sub1 row-index) j)
             (vref row-index (sub1 j))
             (vref (sub1 row-index) (sub1 j))))
      (vset! row-index j (+ cost min-previous-cost))))

  (* earth-radius (vref (sub1 n) (sub1 m))))

;; Same as `dtw`, but uses less memory (and not quadratic to the size of the
;; input), will also be slightly faster, but still runs in quadratic time.
(: dtw/memory-efficient (-> (Vectorof Integer) (Vectorof Integer) Real))
(define (dtw/memory-efficient source target)
  (define-values (n m unit-distance)
    (prepare-source-and-target source target))

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

  ;; Fill the costs in the CURRENT row (at ROW-INDEX), using data from the
  ;; PREVIOUS row.  CURRENT and PREVIOUS are one of the THIS-ROW and
  ;; PREVIOUS-ROW vectors (and they swap around each iteration).
  (: fill-row-with-costs (-> Integer FlVector FlVector Void))
  (define (fill-row-with-costs row-index current previous)
    ;; First element in the current vector only looks "above"
    (define cost (unit-distance row-index 0))
    (define min-previous-cost (flvector-ref previous 0))
    (flvector-set! current 0 (+ cost min-previous-cost))
    (for ([j (in-range 1 m)])
      (define cost (unit-distance row-index j))
      (define min-previous-cost
        (min (flvector-ref previous j)
             (flvector-ref current (sub1 j))
             (flvector-ref previous (sub1 j))))
      (flvector-set! current j (+ cost min-previous-cost))))

  ;; Bootstrap the algorithm -- setup the first element in row 0 to the
  ;; distance between the first geoids in each path and fill in row 1.
  (flvector-set! previous-row 0 (unit-distance 0 0))
  (for ([j (in-range 1 m)])
    (define cost (unit-distance 0 j))
    (define min-previous-cost (flvector-ref previous-row (sub1 j)))
    (flvector-set! previous-row j (+ cost min-previous-cost)))

  ;; Fill in the remaining rows always swapping between the two rows and
  ;; returning the last one that we filled in.
  (define last-line
    (ann (let loop ([row-index : Integer 1]
                    [this-row : FlVector this-row]
                    [previous-row : FlVector previous-row])
           (if (< row-index n)
               (begin
                 (fill-row-with-costs row-index this-row previous-row)
                 (loop (add1 row-index) previous-row this-row))
               ;; Note that we swapped lines already, so we need to return the
               ;; previous one.
               previous-row))
         FlVector))
  (* earth-radius (flvector-ref last-line (sub1 m))))

;; Same as `dtw/memory-efficient`, but uses a sliding window for pairing up
;; costs, instead of matching the entire target vector -- this should produce
;; the same result as the other dtw algorithms, but only if source and target
;; don't deviate too much (the window size) from each other.  This should run
;; faster when source and target are roughly the same length and very large.
;;
(: dtw/window (-> (Vectorof Integer) (Vectorof Integer) Real))
(define (dtw/window source target)
  (define-values (n m unit-distance)
    (prepare-source-and-target source target))
  (define w (get-window-size n m))
  (define nhalf-w (- (exact-round (/ w 2)))) ; negative half window...

  (define this-row (make-flvector w +inf.0))
  (define previous-row (make-flvector w +inf.0))

  ;; Same as fill-row-with-costs from dtw/memory-efficient, but account for
  ;; the fact that the current row (which is now a window) might be offset
  ;; from the previous row.  Also, this version returns the index of the
  ;; element with the smallest cost.
  (: fill-row-with-costs (-> Integer FlVector Integer FlVector Integer Integer))
  (define (fill-row-with-costs row-index current current-offset previous previous-offset)
    (define bro (- current-offset previous-offset)) ; between-row-offset
    (let ([min-previous-cost
           (cond ((= bro 0)
                  (flvector-ref previous 0))
                 ((> bro 0)
                  (min (flvector-ref previous (sub1 bro))
                       (flvector-ref previous bro)))
                 (else ; (< bro 0), should not happen
                  0))]
          [cost (unit-distance row-index (+ 0 current-offset))])
      (flvector-set! current 0 (+ cost min-previous-cost)))

    (for/fold ([min-cost : Flonum (flvector-ref current 0)]
               [min-cost-index : Integer 0]
               #:result min-cost-index)
              ([j (in-range 1 w)])

      (define min-previous-cost
        (let ([j-bro (+ j bro)])
          (cond ((< 0 j-bro w)
                 (min (flvector-ref previous j-bro)
                      (flvector-ref current (sub1 j))
                      (flvector-ref previous (sub1 j-bro))))
                ((< -1 (sub1 j-bro) w)
                 (min (flvector-ref current (sub1 j))
                      (flvector-ref previous (sub1 j-bro))))
                (else
                 (flvector-ref current (sub1 j))))))
      (define cost (unit-distance row-index (+ j current-offset)))
      (define ncost (+ cost min-previous-cost))
      (flvector-set! current j ncost)
      (if (< ncost min-cost)
          (values ncost j)
          (values min-cost min-cost-index))))

  (flvector-set! previous-row 0 (unit-distance 0 0))
  (for ([j (in-range 1 w)])
    (define min-previous-cost (flvector-ref previous-row (sub1 j)))
    (define cost (unit-distance 0 j))
    (flvector-set! previous-row j (+ cost min-previous-cost)))

  (define last-line
    (ann (let loop ([row-index : Integer 1]
                    [this-row : FlVector this-row]
                    [this-row-offset : Integer 0]
                    [previous-row : FlVector previous-row]
                    [previous-row-offset : Integer 0])
           (if (< row-index n)
               (let ([min-cost-index
                      (fill-row-with-costs
                       row-index
                       this-row this-row-offset
                       previous-row previous-row-offset)])
                 (loop
                  (add1 row-index)
                  previous-row
                  ;; Only move the row offset forward, if we move it
                  ;; backwards, we land on the bottom triangle of the row
                  ;; costs, which are always lower and the algorithm will
                  ;; produce unrealistically low costs.
                  (min (max this-row-offset (+ this-row-offset min-cost-index nhalf-w))
                       (- m w))
                  this-row
                  this-row-offset))
               ;; Note that we swapped lines already, so we need to return the
               ;; previous one.
               previous-row))
         FlVector))
  (* earth-radius (flvector-ref last-line (sub1 w))))

(define waypoint-alignment-cost dtw/memory-efficient)
