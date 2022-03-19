#lang typed/racket/base

;; SPDX-License-Identifier: LGPL-3.0-or-later
;; geodesy.rkt -- calculate distances on the earth ellipsoid
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

(require math/flonum
         racket/math)

(provide (all-defined-out))


;;............................................................ utilities ....

;; A small number, used to determine if two numbers are practically equal.
(: ε Flonum)
(define ε (- (flnext 1.0) 1.0))

(: 2π Flonum)
(define 2π (* 2.0 pi))

;; Convert a Number to a Flonum.  If the number is complex, take the real part
;; -- this works for trigonometry, where we might get small imaginary parts at
;; the "edges", which we wish to discard...
(: as-flonum (-> Number Flonum))
(define (as-flonum x)
  (real->double-flonum (if (complex? x) (real-part x) x)))

;; Wrap an angle in radians to the 0 .. 2π range.  This works for all numbers,
;; positive or negative, no matter how large...
(: wrap2π (-> Flonum Flonum))
(define (wrap2π x)
  (+ x (* 2π (ceiling (/ (- x) 2π)))))

;; Wrap an angle in radians to the 0 .. 2π range and shift it by SHIFT amount.
;; If SHIFT is -π, the resulting angle will be in the -π .. π range.
(: wrap2π/shift (-> Flonum Flonum Flonum))
(define (wrap2π/shift x shift)
  (+ x (* 2π (ceiling (/ (- shift x) 2π)))))


;;............................................................ ellipsoid ....

;; Define parameters for an ellipsoid on which distances are calculated -- the
;; Earth is approximated to an ellipsoid, rather than a sphere, as this
;; produces more accurate results.
(struct ellipsoid ([major : Flonum]
                   [minor : Flonum]
                   [flattening : Flonum]
                   [second-eccentricity : Flonum])
  #:transparent)

;; Create an ellipsoid given its semi-major and semi-minor axes.
(: make-ellipsoid (-> Nonnegative-Flonum Nonnegative-Flonum ellipsoid))
(define (make-ellipsoid major minor)
  (ellipsoid major minor 
    (cast (/ (- major minor) major) Flonum)
    (cast (/ (- (sqr major) (sqr minor)) (sqr minor)) Flonum)))

;; The "WGS84" ellipsoid, used by the GPS system...
;; https://en.wikipedia.org/wiki/World_Geodetic_System
(define wgs84 (make-ellipsoid 6378137.0 6356752.314245))


;;............................................................. vincenty ....

;; "first geodetic problem": Find the destination point, given a point
;; (latitude, longitude), a direction (azimuth) and a distance.
;;
;; Arguments: λ - longitude in radians, φ - latitude in radians, α - bearing
;; in radians, 0 being north and turning clockwise, s is distance in meters.
;;
;; Returns the longitude, latitude and final bearing of the destination point.
;; All angles are in radians.
;;
;; See also: https://en.wikipedia.org/wiki/Geodesy,
;;
;; Implementation based on
;; https://movable-type.co.uk/scripts/latlong-vincenty.html

(: vincenty-direct (->* (Real Real Real Real)
                        (#:ellipsoid ellipsoid
                         #:max-iterations Positive-Integer
                         #:precision Positive-Real)
                        (Values Real Real Real)))
(define (vincenty-direct λ1 φ1 α1 s
                         #:ellipsoid (e wgs84)
                         #:max-iterations (max-iterations 1000)
                         #:precision (precision 1e-12))
  (define a (ellipsoid-major e))
  (define b (ellipsoid-minor e))
  (define f (ellipsoid-flattening e))
  (define 2e (ellipsoid-second-eccentricity e))

  (define sinα1 (sin α1))
  (define cosα1 (cos α1))
  (define tanU1 (as-flonum (* (- 1.0 f) (tan φ1))))
  (define cosU1 (/ 1.0 (sqrt (+ 1.0 (sqr tanU1)))))
  (define sinU1 (* tanU1 cosU1))
  ;; σ1 = angular distance on the sphere from the equator to P1
  (define σ1 (atan tanU1 cosα1))
  (define sinα (as-flonum (* cosU1 sinα1))) ; α - azimuth of the geodesic at the equator
  (define cos²α (- 1.0 (sqr sinα)))
  (define u² (* cos²α 2e))

  (define A
    (+ 1.0
       (let loop : Flonum
            ([result : Flonum 0.0]
             [x : Flonum u²]
             [coefficients : (Listof Flonum) '(4096.0 -768.0 320.0 -175.0)])
         (if (null? coefficients)
             (/ result 16384.0)
             (loop (+ result (* (car coefficients) x))
                   (* x u²)
                   (cdr coefficients))))))
  (define B
    (let loop : Flonum
         ([result : Flonum 0.0]
          [x : Flonum u²]
          [coefficients : (Listof Flonum) '(256.0 -128.0 74.0 -47.0)])
      (if (null? coefficients)
          (/ result 1024.0)
          (loop (+ result (* (car coefficients) x))
                (* x u²)
                (cdr coefficients)))))

  (: σ Flonum)
  (: sinσ Flonum)
  (: cosσ Flonum)
  (: cos2σₘ Flonum)
  (define-values (σ sinσ cosσ cos2σₘ)
    (let loop ([σ : Flonum (as-flonum (/ s (* b A)))]
               [remaining-iterations : Integer max-iterations])
      (define cos2σₘ (cos (+ (* 2.0 σ1) σ)))
      (define sinσ (sin σ))
      (define cosσ (cos σ))
      (define Δσ
        (let ([cos²2σₘ (sqr cos2σₘ)]
              [sin²σ (sqr sinσ)])
          (* B
             sinσ
             (+ cos2σₘ
                (* (/ B 4.0)
                   (- (* cosσ (+ -1.0 (* 2.0 cos²2σₘ)))
                      (* (/ B 6.0) cos2σₘ (+ -3.0 (* 4.0  sin²σ)) (+ -3.0  (* 4.0 cos²2σₘ)))))))))
      (define new-σ (+ (/ s (* b A)) Δσ))
      (cond ((< (abs (- new-σ σ)) precision)
             (values σ sinσ cosσ cos2σₘ))
            ((zero? remaining-iterations)
             (error "vincenty-direct: formula failed to converge"))
            (#t
             (loop new-σ (sub1 remaining-iterations))))))

  (define x (as-flonum (- (* sinU1 sinσ) (* cosU1 cosσ cosα1))))
  (define φ2 (atan (as-flonum (+ (* sinU1 cosσ) (* cosU1 sinσ cosα1)))
                   (as-flonum (* (- 1.0 f) (sqrt (+ (sqr sinα) (sqr x)))))))
  (define λ (atan (* sinσ sinα1) (as-flonum (- (* cosU1 cosσ) (* sinU1 sinσ cosα1)))))

  (define C (* (/ f 16.0) cos²α (+ 4.0 (* f (- 4.0 (* 3 cos²α))))))
  (define L (- λ (* (- 1.0 C)
                    f
                    sinα
                    (+ σ (* C sinσ (+ cos2σₘ (* C cosσ (+ -1.0 (* 2.0 (sqr cos2σₘ))))))))))
  (define λ2 (+ λ1 L))
  (define α2 (atan sinα (- x)))

  (values (wrap2π/shift λ2 (- pi)) φ2 (wrap2π α2)))

;; The second geodetic problem: given two points (latitude, longitude), find
;; the distance, initial and final bearings for the great circle path between
;; them.
;;
;; λ - longitude in radians, φ - latitude in radians.
;;
;; Returns 3 values: distance, initial bearing and final bearing -- returned
;; angles are in radians 0 being north and moving clockwise.
;; See also: https://en.wikipedia.org/wiki/Geodesy,
;;
;; Implementation based on
;; https://movable-type.co.uk/scripts/latlong-vincenty.html

(: vincenty-inverse (->* (Real Real Real Real)
                         (#:ellipsoid ellipsoid
                          #:max-iterations Positive-Integer
                          #:precision Positive-Real)
                         (Values Real Real Real)))
(define (vincenty-inverse λ1 φ1 λ2 φ2
                          #:ellipsoid (e wgs84)
                          #:max-iterations (max-iterations 1000)
                          #:precision (precision 1e-12))

  (define a (ellipsoid-major e))
  (define b (ellipsoid-minor e))
  (define f (ellipsoid-flattening e))
  (define 2e (ellipsoid-second-eccentricity e))

  (define L (as-flonum (- λ2 λ1)))

  (define tan-U1 (as-flonum (* (- 1.0 f) (tan φ1))))
  (define cos-U1 (as-flonum (/ 1.0 (sqrt (+ 1.0 (* tan-U1 tan-U1))))))
  (define sin-U1 (as-flonum (* tan-U1 cos-U1)))

  (define tan-U2 (as-flonum (* (- 1.0 f) (tan φ2))))
  (define cos-U2 (as-flonum (/ 1 (sqrt (+ 1 (* tan-U2 tan-U2))))))
  (define sin-U2 (as-flonum (* tan-U2 cos-U2)))

  (define antipodal? (or (> (abs L) (/ pi 2))
                         (> (abs (- φ2 φ1)) (/ pi 2))))

  (: cos²α Flonum)
  (: sin²σ Flonum)
  (: cos2σₘ Flonum)
  (: sinσ Flonum)
  (: cosσ Flonum)
  (: σ Flonum)
  (: sinλ Flonum)
  (: cosλ Flonum)
  (define-values
    (sinλ cosλ σ sinσ cosσ sin²σ cos2σₘ cos²α)
    (let loop ([λ : Flonum L]
               [sinλ : Flonum (as-flonum (sin L))]
               [cosλ : Flonum (as-flonum (cos L))]
               [σ : Flonum (if antipodal? pi 0.0)]
               [sinσ : Flonum 0.0]
               [cosσ : Flonum (if antipodal? -1.0 1.0)]
               [cos2σₘ : Flonum 1.0]
               [cos²α : Flonum 1.0]
               [remaining-iterations : Integer max-iterations])

      (define sin²σ (+ (sqr (* cos-U2 sinλ))
                       (sqr (- (* cos-U1 sin-U2) (* sin-U1 cos-U2 cosλ)))))
      (if (< sin²σ ε)
          (begin
            ;; exit early, co-incident, antipodal points
            (values sinλ cosλ σ sinσ cosσ 0.0 cos2σₘ cos²α))
          (let* ([sinσ : Flonum (as-flonum (sqrt sin²σ))]
                 [cosσ : Flonum (+ (* sin-U1 sin-U2) (* cos-U1 cos-U2 cosλ))]
                 [σ : Flonum (as-flonum (atan sinσ cosσ))]
                 [sinα : Flonum (/ (* cos-U1 cos-U2 sinλ) sinσ)]
                 [cos²α : Flonum (- 1.0 (* sinα sinα))]
                 [cos2σₘ : Flonum (if (zero? cos²α)
                                      0.0
                                      (- cosσ (/ (* 2.0 sin-U1 sin-U2) cos²α)))]
                 [C : Flonum (* (/ f 16.0) cos²α (+ 4 (* f (- 4 (* 3 cos²α)))))]
                 [new-λ : Flonum
                        (+ L
                           (* (- 1.0 C)
                              f
                              sinα
                              (+ σ
                                 (* C sinσ (+ cos2σₘ (* C cosσ (+ -1.0 (* 2.0 cos2σₘ cos2σₘ))))))))])

            (define check (if antipodal? (- (abs new-λ) pi) new-λ))
            (when (> check pi)
              (error (format "new-λ > pi, new-λ : ~a, antipodal? ~a" new-λ antipodal?)))

            (if (and (> remaining-iterations 0)
                     (> (abs (- new-λ λ)) precision))
                (loop new-λ (sin new-λ) (cos new-λ) σ sinσ cosσ cos2σₘ cos²α (sub1 remaining-iterations))
                (values sinλ cosλ σ sinσ cosσ sin²σ cos2σₘ cos²α))))))

  (define u² (* cos²α 2e))

  (define A
    (+ 1.0
       (let loop : Flonum
            ([result : Flonum 0.0]
             [x : Flonum u²]
             [coefficients : (Listof Flonum) '(4096.0 -768.0 320.0 -175.0)])
         (if (null? coefficients)
             (/ result 16384.0)
             (loop (+ result (* (car coefficients) x))
                   (* x u²)
                   (cdr coefficients))))))
  (define B
    (let loop : Flonum
         ([result : Flonum 0.0]
          [x : Flonum u²]
          [coefficients : (Listof Flonum) '(256.0 -128.0 74.0 -47.0)])
      (if (null? coefficients)
          (/ result 1024.0)
          (loop (+ result (* (car coefficients) x))
                (* x u²)
                (cdr coefficients)))))

  (define Δσ
    (let ([cos²2σₘ (* cos2σₘ cos2σₘ)])
      (* B
         sinσ
         (+ cos2σₘ
            (* (/ B 4)
               (- (* cosσ (+ -1.0 (* 2 cos²2σₘ)))
                  (* (/ B 6) cos2σₘ (+ -3.0 (* 4.0  sin²σ)) (+ -3.0  (* 4 cos²2σₘ)))))))))

  (define s (* b A (- σ Δσ)))           ; length of the geodesic
  (define α1 (if (< (abs sin²σ) ε)
                 0.0
                 (atan (* cos-U2 sinλ) (- (* cos-U1 sin-U2) (* sin-U1 cos-U2 cosλ)))))
  (define α2 (if (< (abs sin²σ) ε)
                 pi
                 (atan (* cos-U1 sinλ) (- (* cos-U1 sin-U2 cosλ) (* sin-U1 cos-U2)))))

  (if (< s ε)
      (values s +nan.0 +nan.0)
      (values s (wrap2π α1) (wrap2π α2))))


;;............................................................... sphere ....

;; Distances and destination points, assuming the earth is a sphere.
;; Implementation based on routines from:
;; https://movable-type.co.uk/scripts/latlong.html

;; Radius from https://en.wikipedia.org/wiki/Reference_ellipsoid
(: earth-radius Real)
(define earth-radius 6371088.0)

(: spherical-distance (-> Real Real Real Real Real))
(define (spherical-distance λ1 φ1 λ2 φ2)
  (: Δφ/2 Flonum)
  (define Δφ/2 (as-flonum (/ (- φ2 φ1) 2.0)))
  (: Δλ/2 Flonum)
  (define Δλ/2 (as-flonum (/ (- λ2 λ1) 2.0)))
  (define a (+ (sqr (sin Δφ/2))
               (* (cos φ1) (cos φ2) (sqr (as-flonum (sin Δλ/2))))))
  (define c (* 2.0 (atan (as-flonum (sqrt a))
                         (as-flonum (sqrt (- 1 a))))))
  (* earth-radius c))

(: spherical-bearing (-> Real Real Real Real Real))
(define (spherical-bearing λ1 φ1 λ2 φ2)
  (define Δλ (- λ2 λ1))
  (define y (* (sin Δλ) (cos φ2)))
  (define x (- (* (cos φ1) (sin φ2))
               (* (sin φ1) (cos φ2) (cos Δλ))))
  (wrap2π (as-flonum (atan y x))))

(: spherical-destination (-> Real Real Real Real (values Real Real)))
(define (spherical-destination λ φ bearing distance)
  (define δ (/ distance earth-radius))
  (define φ-dest
    (asin (+ (* (sin φ) (cos δ))
             (* (cos φ) (sin δ) (cos bearing)))))
  (define λ-dest
    (+ λ
       (atan (* (sin bearing) (sin δ) (cos φ))
             (- (cos δ) (* (sin φ) (sin φ-dest))))))
  (values (wrap2π/shift (as-flonum λ-dest) (- pi)) φ-dest))

(: spherical-intermediate-point (-> Real Real Real Real Real (values Real Real)))
(define (spherical-intermediate-point λ1 φ1 λ2 φ2 f)
  (define distance (spherical-distance λ1 φ1 λ2 φ2))
  (define δ (/ distance earth-radius))
  (define a (/ (sin (* (- 1 f) δ)) (sin δ)))
  (define b (/ (sin (* f δ)) (sin δ)))
  (define x (as-flonum (+ (* a (cos φ1) (cos λ1))
                          (* b (cos φ2) (cos λ2)))))
  (define y (as-flonum (+ (* a (cos φ1) (sin λ1))
                          (* b (cos φ2) (sin λ2)))))
  (define z (+ (* a (sin φ1)) (* b (sin φ2))))
  (values (wrap2π/shift (atan y x) (- pi))
          (atan z (as-flonum (sqrt (+ (sqr x) (sqr y)))))))

(: spherical-midway-point (-> Real Real Real Real (values Real Real)))
(define (spherical-midway-point λ1 φ1 λ2 φ2)
  (define Δλ (- λ2 λ1))
  (define Bx (as-flonum (* (cos φ2) (cos Δλ))))
  (define By (as-flonum (* (cos φ2) (sin Δλ))))
  (values (wrap2π/shift (as-flonum (+ λ1 (atan By (+ (cos φ1) Bx)))) (- pi))
          (as-flonum (atan (+ (sin φ1) (sin φ2))
                           (as-flonum (sqrt (+ (sqr (+ (cos φ1) Bx))
                                               (sqr By))))))))


;;.................................................................. API ....

;; Convert input coordinates from angles to radians and re-order them.
;; Internal functions use the "longitude, latitude" argument order, but API
;; use "latitude, longitude"
(: ->radians (-> Real Real  (U 'degrees 'radians)
                 (Values Flonum Flonum)))
(define (->radians lat lon angle-mode)
  (if (eq? angle-mode 'degrees)
      (values
       (as-flonum (degrees->radians lon))
       (as-flonum (degrees->radians lat)))
      (values
       (as-flonum lon) (as-flonum lat))))

;; Angle mode (degrees or radians) for input parameters of the API functions.
(: geodesy-angle-mode (Parameter (U 'degrees 'radians)))
(define geodesy-angle-mode (make-parameter 'degrees))

;; Ellipsoid to use for the calculations.  #f means use a sphere
(: geodesy-ellipsoid (Parameter (U #f ellipsoid)))
(define geodesy-ellipsoid (make-parameter wgs84))


(: distance-between (->* (Real Real Real Real)
                         (#:angle-mode (U 'degrees 'radians)
                          #:ellipsoid (U #f ellipsoid))
                         Real))
(define (distance-between lat1 lon1 lat2 lon2
                          #:angle-mode (mode (geodesy-angle-mode))
                          #:ellipsoid (ellipsoid (geodesy-ellipsoid)))
  (define-values (λ1 φ1) (->radians lat1 lon1 mode))
  (define-values (λ2 φ2) (->radians lat2 lon2 mode))
  (if ellipsoid
      (let-values ([(distance _initial-bearing _final-bearing)
                    (vincenty-inverse λ1 φ1 λ2 φ2 #:ellipsoid ellipsoid)])
        distance)
      (spherical-distance λ1 φ1 λ2 φ2)))

(: initial-bearing (->* (Real Real Real Real)
                        (#:angle-mode (U 'degrees 'radians)
                         #:ellipsoid (U #f ellipsoid))
                        Real))
(define (initial-bearing lat1 lon1 lat2 lon2
                         #:angle-mode (mode (geodesy-angle-mode))
                         #:ellipsoid (ellipsoid (geodesy-ellipsoid)))
  (define-values (λ1 φ1) (->radians lat1 lon1 mode))
  (define-values (λ2 φ2) (->radians lat2 lon2 mode))
  (define b
    (if ellipsoid
        (let-values ([(_distance initial-bearing _final-bearing)
                      (vincenty-inverse λ1 φ1 λ2 φ2 #:ellipsoid ellipsoid)])
          initial-bearing)
        (spherical-bearing λ1 φ1 λ2 φ2)))
  (if (eq? mode 'degrees) (radians->degrees b) b))

(: final-bearing (->* (Real Real Real Real)
                      (#:angle-mode (U 'degrees 'radians)
                       #:ellipsoid (U #f ellipsoid))
                      Real))
(define (final-bearing lat1 lon1 lat2 lon2
                       #:angle-mode (mode (geodesy-angle-mode))
                       #:ellipsoid (ellipsoid (geodesy-ellipsoid)))
  (define-values (λ1 φ1) (->radians lat1 lon1 mode))
  (define-values (λ2 φ2) (->radians lat2 lon2 mode))
  (define b
    (if ellipsoid
        (let-values ([(_distance _initial-bearing final-bearing)
                      (vincenty-inverse λ1 φ1 λ2 φ2 #:ellipsoid ellipsoid)])
          final-bearing)
        (wrap2π (- (spherical-bearing λ2 φ2 λ1 φ1) pi))))
  (if (eq? mode 'degrees) (radians->degrees b) b))

(: destination-point (->* (Real Real Real Real)
                          (#:angle-mode (U 'degrees 'radians)
                           #:ellipsoid (U #f ellipsoid))
                          (values Real Real)))
(define (destination-point lat lon bearing distance
                           #:angle-mode (mode (geodesy-angle-mode))
                           #:ellipsoid (ellipsoid (geodesy-ellipsoid)))
  (define-values (λ φ) (->radians lat lon mode))
  (define α (if (eq? mode 'degrees) (degrees->radians bearing) bearing))
  (if ellipsoid
      (let-values ([(λ-dest φ-dest _final-bearing)
                    (vincenty-direct λ φ α distance #:ellipsoid ellipsoid)])
        (if (eq? mode 'degrees)
            (values (radians->degrees φ-dest) (radians->degrees λ-dest))
            (values φ-dest λ-dest)))
      (let-values ([(λ-dest φ-dest) (spherical-destination λ φ α distance)])
        (if (eq? mode 'degrees)
            (values (radians->degrees φ-dest) (radians->degrees λ-dest))
            (values φ-dest λ-dest)))))

(: midway-point (->* (Real Real Real Real)
                     (#:angle-mode (U 'degrees 'radians)
                      #:ellipsoid (U #f ellipsoid))
                     (values Real Real)))
(define (midway-point lat1 lon1 lat2 lon2
                      #:angle-mode (mode (geodesy-angle-mode))
                      #:ellipsoid (ellipsoid (geodesy-ellipsoid)))
  (define-values (λ1 φ1) (->radians lat1 lon1 mode))
  (define-values (λ2 φ2) (->radians lat2 lon2 mode))
  (define-values (λm φm)
    (if ellipsoid
        (let-values ([(distance initial-bearing _final-bearing)
                      (vincenty-inverse λ1 φ1 λ2 φ2 #:ellipsoid ellipsoid)])
          (define-values (λm φm _intermediate_bearig)
            (vincenty-direct λ1 φ1 initial-bearing (* 0.5 distance) #:ellipsoid ellipsoid))
          (values λm φm))
        (spherical-midway-point λ1 φ1 λ2 φ2)))
  (if (eq? mode 'degrees)
      (values (radians->degrees φm) (radians->degrees λm))
      (values φm λm)))
