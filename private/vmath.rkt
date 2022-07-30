#lang typed/racket/base

;; SPDX-License-Identifier: LGPL-3.0-or-later
;;
;; vmath.rkt -- 3D vector operations
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

(require math/flonum math/bigfloat)
(provide (all-defined-out))

;; 3D Flonum vectors are the basis for the tiling operations (segment
;; intersections).  Here we implement the usual operations on them...

(define-type Vector3 FlVector)

;; A "greater than" operation to compare and order 3D Flonum vectors.  Used to
;; order the triangle edges for `symbolically-perturbed-sign`
(: v3gt (-> Vector3 Vector3 Boolean))
(define (v3gt a b)
  (define a0 (flvector-ref a 0))
  (define a1 (flvector-ref a 1))
  (define a2 (flvector-ref a 2))
  (define b0 (flvector-ref b 0))
  (define b1 (flvector-ref b 1))
  (define b2 (flvector-ref b 2))
  (cond ((> a0 b0) #t)
        ((< a0 b0) #f)
        ((> a1 b1) #t)                  ; a0 == b0
        ((< a1 b1) #f)
        ((> a2 b2) #t)                  ; a1 == b1
        ((< a2 b2) #f)
        ;; The two vectors are equal, so they are not greater than each
        ;; other...
        (#t #f)))

(: v3+ (case->
        [Vector3 -> Vector3]
        [Vector3 Vector3 -> Vector3]
        [->* (Vector3 Vector3) #:rest Vector3 Vector3]))
(define v3+
  (case-lambda
    [(a) a]
    [(a b) (flvector+ a b)]
    [(a b . rest)
     (let loop : Vector3
          ([result : Vector3 (flvector+ a b)]
           [remaining : (Listof Vector3) rest])
       (if (null? remaining)
           result
           (loop (flvector+ result (car remaining)) (cdr remaining))))]))

(: v3- (case->
        [Vector3 -> Vector3]
        [Vector3 Vector3 -> Vector3]
        [->* (Vector3 Vector3) #:rest Vector3 Vector3]))
(define v3-
  (case-lambda
    [(a) (flvector- a)]
    [(a b) (flvector- a b)]
    [(a b . rest)
     (let loop : Vector3
          ([result : Vector3 (flvector- a b)]
           [remaining : (Listof Vector3) rest])
       (if (null? remaining)
           result
           (loop (flvector- result (car remaining)) (cdr remaining))))]))

(: v3scale (-> Vector3 Flonum Vector3))
(define (v3scale a s)
  (flvector-scale a s))

(: v3unit (-> Vector3 Vector3))
(define (v3unit a)
  (flvector-scale a (/ 1.0 (v3len a))))

(: v3len² (-> Vector3 Flonum))
(define (v3len² a)
  (dot-product a a))

(: v3len (-> Vector3 Flonum))
(define (v3len a)
  (assert (sqrt (v3len² a)) flonum?))

(: dot-product (-> Vector3 Vector3 Flonum))
(define (dot-product a b)
  (+
   (* (flvector-ref a 0) (flvector-ref b 0))
   (* (flvector-ref a 1) (flvector-ref b 1))
   (* (flvector-ref a 2) (flvector-ref b 2))))


(: cross-product (-> Vector3 Vector3 Vector3))
(define (cross-product a b)
  (flvector (- (* (flvector-ref a 1) (flvector-ref b 2))
               (* (flvector-ref a 2) (flvector-ref b 1)))
            (- (* (flvector-ref a 2) (flvector-ref b 0))
               (* (flvector-ref a 0) (flvector-ref b 2)))
            (- (* (flvector-ref a 0) (flvector-ref b 1))
               (* (flvector-ref a 1) (flvector-ref b 0)))))

;; https://github.com/google/s2geometry/blob/0c4c460bdfe696da303641771f9def900b3e440f/src/s2/s2pointutil.cc#L57
;;
;; This version of cross-product works when the vectors a and b are close to
;; each other, that is, (a + b) or (a - b) are close to zero.  I don't know
;; how or why it works, but it does.  The implementation is inspired from the
;; link above.
(: robust-cross-product (-> Vector3 Vector3 Vector3))
(define (robust-cross-product a b)
  (define ua (v3unit a))
  (define ub (v3unit b))
  (cross-product (v3+ ub ua) (v3- ub ua)))


;; 3D Bigfloat vectors are used for triangle operations for such thin
;; triangles that normal flonums cannot handle them.  Here we define the usual
;; 3D vector operations for them.

(define-type BfVector3 (Vector Bigfloat Bigfloat Bigfloat))

(: v3->bfv3 (-> Vector3 BfVector3))
(define (v3->bfv3 v)
  (vector (bf (flvector-ref v 0))
          (bf (flvector-ref v 1))
          (bf (flvector-ref v 2))))

(: bfv3+ (case->
          [BfVector3 -> BfVector3]
          [BfVector3 BfVector3 -> BfVector3]
          [->* (BfVector3 BfVector3) #:rest BfVector3 BfVector3]))
(define bfv3+
  (case-lambda
    [(a) a]
    [(a b)
     (vector (bf+ (vector-ref a 0) (vector-ref b 0))
             (bf+ (vector-ref a 1) (vector-ref b 1))
             (bf+ (vector-ref a 2) (vector-ref b 2)))]
    [(a b . rest)
     (let loop : BfVector3
          ([result : BfVector3 (bfv3+ a b)]
           [remaining : (Listof BfVector3) rest])
       (if (null? remaining)
           result
           (loop (bfv3+ result (car remaining)) (cdr remaining))))]))


(: bfv3- (case->
          [BfVector3 -> BfVector3]
          [BfVector3 BfVector3 -> BfVector3]
          [->* (BfVector3 BfVector3) #:rest BfVector3 BfVector3]))
(define bfv3-
  (case-lambda
    [(a) (vector (bf- (vector-ref a 0))
                 (bf- (vector-ref a 1))
                 (bf- (vector-ref a 2)))]
    [(a b) (vector (bf- (vector-ref a 0) (vector-ref b 0))
                   (bf- (vector-ref a 1) (vector-ref b 1))
                   (bf- (vector-ref a 2) (vector-ref b 2)))]
    [(a b . rest)
     (let loop : BfVector3
          ([result : BfVector3 (bfv3- a b)]
           [remaining : (Listof BfVector3) rest])
       (if (null? remaining)
           result
           (loop (bfv3- result (car remaining)) (cdr remaining))))]))

(: bfv3scale (-> BfVector3 Bigfloat BfVector3))
(define (bfv3scale a s)
  (vector (bf* (vector-ref a 0) s)
          (bf* (vector-ref a 1) s)
          (bf* (vector-ref a 2) s)))

(: bfv3unit (-> BfVector3 BfVector3))
(define (bfv3unit a)
  (bfv3scale a (bf/ 1.bf (bfv3len a))))

(: bfv3len² (-> BfVector3 Bigfloat))
(define (bfv3len² a)
  (bf-dot-product a a))

(: bfv3len (-> BfVector3 Bigfloat))
(define (bfv3len a)
  (bfsqrt (bfv3len² a)))

(: bf-dot-product (-> BfVector3 BfVector3 Bigfloat))
(define (bf-dot-product a b)
  (bf+
   (bf* (vector-ref a 0) (vector-ref b 0))
   (bf* (vector-ref a 1) (vector-ref b 1))
   (bf* (vector-ref a 2) (vector-ref b 2))))

(: bf-cross-product (-> BfVector3 BfVector3 BfVector3))
(define (bf-cross-product a b)
  (vector (bf- (bf* (vector-ref a 1) (vector-ref b 2))
               (bf* (vector-ref a 2) (vector-ref b 1)))
          (bf- (bf* (vector-ref a 2) (vector-ref b 0))
               (bf* (vector-ref a 0) (vector-ref b 2)))
          (bf- (bf* (vector-ref a 0) (vector-ref b 1))
               (bf* (vector-ref a 1) (vector-ref b 0)))))

;; Triangle sign routines are implemented based on the Google S2 library
;; implementations (see the equivalent functions in s2predicates.cc, for
;; details of how and why they work).
;;
;; The entry routines are `sign` and `sign/c` (the latter just uses a
;; cross-product, if it is already available).  These routines will use more
;; specific routines, such as `triage-sign`, `stable-sign`, `exact-sign` and
;; `symbolically-perturbed-sign`, which are slower and slower but can have
;; better precision.  `symbolically-perturbed-sign` is a special one, since it
;; will determinisitally assign a sign to *every* triangle even if the A, B, C
;; edges are co-linear.  By the way, sign determines the winding order of a
;; triangle determined by the vertices A, B and C and returns -1 for clockwise
;; winding, 1 for counter-clockwise and 0 if two vertices are identical.

;; Maximum error for the determinant of the triangle, as determined by
;; `triage-sign` -- for determinants smaller than this, triage sign returns 0,
;; meaning it cannot determine the winding order.
(define max-triage-sign-error (* 1.8274 epsilon.0))

(: triage-sign (-> Vector3 Vector3 Vector3 Vector3 Integer))
(define (triage-sign a b c a-cross-b)
  (define det (dot-product a-cross-b c))
  (cond ((> det max-triage-sign-error) 1)
        ((< det (- max-triage-sign-error)) -1)
        (#t 0)))

;; Errors for the `stable-sign` determinant determinantion, for numbers
;; smaller than these, `stable-sign` returns 0, meaning it cannot determine
;; the winding order.
(define stable-sign-error-multipler (* 3.2321 epsilon.0))
(define stable-sign-no-underflow-error (* stable-sign-error-multipler +min.0))

(: stable-sign (-> Vector3 Vector3 Vector3 Integer))
(define (stable-sign a b c)
  (define ab (v3- b a))
  (define bc (v3- c b))
  (define ca (v3- a c))
  (define ab2 (v3len² ab))
  (define bc2 (v3len² bc))
  (define ca2 (v3len² ca))

  (define-values (det max-error2)
    (cond
      ((and (>= ab2 bc2) (>= ab2 ca2))
       (values (- (dot-product (cross-product ca bc) c))
               (* ca2 bc2)))
      ((>= bc2 ca2)
       (values (- (dot-product (cross-product ab ca) a))
               (* ab2 ca2)))
      (#t
       (values (- (dot-product (cross-product bc ab) b))
               (* bc2 ab2)))))
  (define max-error
    (* stable-sign-error-multipler (assert (sqrt max-error2) flonum?)))

  (cond ((< max-error stable-sign-no-underflow-error)
         0)
        ((<= (abs det) max-error)
         0)
        ((> det 0)
         1)
        (#t
         -1)))

;; Determine the winding order of 3 vertices A, B and C, which are known to be
;; on the same line.  Also A, B, C should be ordered according to `v3gt` -- in
;; fact, this function should only be called from `exact-sign`.
(: symbolically-perturbed-sign (-> BfVector3 BfVector3 BfVector3 BfVector3 Integer))
(define (symbolically-perturbed-sign a b c b-cross-c)

  (: sgn1 (-> Bigfloat Integer))
  (define (sgn1 b) (bigfloat->integer (bfsgn b)))

  (: sgn2 (-> BfVector3 Integer Integer))
  (define (sgn2 v index)
    (sgn1 (vector-ref v index)))

  (let/ec return : Integer
    (let ([sign (sgn2 b-cross-c 2)])
      (unless (zero? sign) (return sign)))
    (let ([sign (sgn2 b-cross-c 1)])
      (unless (zero? sign) (return sign)))
    (let ([sign (sgn2 b-cross-c 0)])
      (unless (zero? sign) (return sign)))

    (let ([sign (sgn1 (bf- (bf* (vector-ref c 0)
                                (vector-ref a 1))
                           (bf* (vector-ref c 1)
                                (vector-ref a 0))))])
      (unless (zero? sign) (return sign)))

    (let ([sign (sgn2 c 0)])
      (unless (zero? sign) (return sign)))

    (let ([sign (- (sgn2 c 1))])
      (unless (zero? sign) (return sign)))

    (let ([sign (sgn1 (bf- (bf* (vector-ref c 2)
                                (vector-ref a 0))
                           (bf* (vector-ref c 0)
                                (vector-ref a 2))))])
      (unless (zero? sign) (return sign)))

    (let ([sign (sgn2 c 2)])
      (unless (zero? sign) (return sign)))

    (let ([sign (sgn1 (bf- (bf* (vector-ref a 0)
                                (vector-ref b 1))
                           (bf* (vector-ref a 1)
                                (vector-ref b 0))))])
      (unless (zero? sign) (return sign)))

    (let ([sign (- (sgn2 b 0))])
      (unless (zero? sign) (return sign)))

    (let ([sign (sgn2 b 1)])
      (unless (zero? sign) (return sign)))

    (let ([sign (sgn2 a 0)])
      (unless (zero? sign) (return sign)))

    1))

;; Determine the winding order of the triangle A, B, C without any errors
;; (uses Bigfloat for the calculations, which is expensive).  If A, B, C are
;; on the same line, calls `symbolically-perturbed-sign` to resolve the
;; winding order.
(: exact-sign (-> Vector3 Vector3 Vector3 Integer))
(define (exact-sign a b c)
  (define perm-sign 1)

  ;; We order the vertices, since symbolically-perturbed-sign expects them
  ;; ordered.  Permuting the vertices changes the sign, and we keep track of
  ;; that in `perm-sign`
  (when (v3gt a b)
    (let ([tmp a])
      (set! a b)
      (set! b tmp)
      (set! perm-sign (- perm-sign))))
  (when (v3gt b c)
    (let ([tmp b])
      (set! b c)
      (set! c tmp)
      (set! perm-sign (- perm-sign))))
  (when (v3gt a b)
    (let ([tmp a])
      (set! a b)
      (set! b tmp)
      (set! perm-sign (- perm-sign))))

  (define xa (v3->bfv3 a))
  (define xb (v3->bfv3 b))
  (define xc (v3->bfv3 c))

  (define xb-cross-xc (bf-cross-product xb xc))
  (define det (bf-dot-product xa xb-cross-xc))

  (define det-sign (bigfloat->integer (bfsgn det)))
  (* perm-sign
     (if (zero? det-sign)
         (symbolically-perturbed-sign xa xb xc xb-cross-xc)
         det-sign)))

;; Return the winding order of the triangle determined by the vertices A, B
;; and C.  Returns 0 ONLY if two of the vertices are the same, otherwise
;; returns -1 for CW ordering and 1 for CCW ordering.  This is a convenience
;; function for when the user has the cross product between A and B already
;; calculated.
(: sign/c (-> Vector3 Vector3 Vector3 Vector3 Integer))
(define (sign/c a b c a-cross-b)
  ;; Return 0 if any two points are the same
  (if (or (equal? a b) (equal? a c) (equal? b c))
      0
      (let ([s0 (triage-sign a b c a-cross-b)])
        (if (zero? s0)
            (let ([s1 (stable-sign a b c)])
              (if (zero? s1)
                  (exact-sign a b c)
                  s1))
            s0))))

;; Return the winding order of the triangle determined by the vertices A, B
;; and C.  Returns 0 ONLY if two of the vertices are the same, otherwise
;; returns -1 for CW ordering and 1 for CCW ordering.
(: sign (-> Vector3 Vector3 Vector3 Integer))
(define (sign a b c)
  (sign/c a b c (cross-product a b)))

;; Return the index of the larges absolute value in the vector V
(: v3-largest-abs-component (-> Vector3 Integer))
(define (v3-largest-abs-component v)
  (define a (abs (flvector-ref v 0)))
  (define b (abs (flvector-ref v 1)))
  (define c (abs (flvector-ref v 2)))
  (if (> a b)
      (if (> a c) 0 2)
      (if (> b c) 1 2)))

;; Return #t if the points A, B, C are ordered counter-clockwise with respect
;; to the point O.
(: ordered-ccw (-> Vector3 Vector3 Vector3 Vector3 Boolean))
(define (ordered-ccw a b c o)
  (define sum
    (+ (if (>= (sign b o a) 0) 1 0)
       (if (>= (sign c o b) 0) 1 0)
       ;; NOTE: last one is '>', NOT '>='
       (if (> (sign a o c) 0) 1 0)))
  (>= sum 2))
