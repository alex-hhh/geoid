#lang racket

;; geoid-tests.rkt -- tests for the geoid package
;;
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


(require rackunit "geoid.rkt" math/statistics math/flonum)

(define earth-radius (->fl 6371000))    ; meters

(define (haversin Θ)
  (fl/ (fl- 1.0 (flcos Θ)) 2.0))

(define (inv-haversin h)
  (fl* 2.0 (flasin (flsqrt h))))

;; Calculate the distance in meters between two map coordinates
(define (map-distance/radians lat1 lon1 lat2 lon2)
  (let ((Δ-lat (fl- lat2 lat1))
        (Δ-lon (fl- lon2 lon1)))
    (let* ((a (fl+ (haversin Δ-lat)
                   (fl* (fl* (flcos lat1) (flcos lat2))
                        (haversin Δ-lon))))
           (c (inv-haversin a)))
      (fl* c earth-radius))))

(define (map-distance/degrees lat1 lon1 lat2 lon2)
  (map-distance/radians
   (degrees->radians lat1)
   (degrees->radians lon1)
   (degrees->radians lat2)
   (degrees->radians lon2)))

(define (pp-statistics stats headline)
  (printf "~%*** ~a~%    sample count: ~a~%    mean:  ~a~%    stddev ~a~%    min    ~a~%    max    ~a~%"
          headline
          (statistics-count stats)
          (statistics-mean stats)
          (statistics-stddev stats)
          (statistics-min stats)
          (statistics-max stats)))

(define hilbert-distance-test-suite
  (test-suite
   "Hilbert Distance"
   (test-case
       "Basic Encoding + Decoding"

     ;; These are all hilbert distances for recurse level 1
     (define hencodings-level1
       '((0 0 0) (1 0 1) (2 0 14) (3 0 15) (0 1 3) (1 1 2)
         (2 1 13) (3 1 12) (0 2 4) (1 2 7) (2 2 8) (3 2 11)
         (0 3 5) (1 3 6) (2 3 9) (3 3 10)))

     (for ([data (in-list hencodings-level1)])
       (match-define (list x y d) data)
       (define e (xy->hilbert-distance 1 x y))
       (check-equal? d e (format "Failed encoding for x = ~a, y = ~a" x y))
       (define-values (dx dy) (hilbert-distance->xy 1 d))
       (check-equal? dx x (format "Failed decoding for ~a" d))
       (check-equal? dy y (format "Failed decoding for ~a" d))))

   (test-case
       "Higher Level Encode + Decode"

     ;; For higher level encodings, we just make sure that the code encodes
     ;; and decodes to the same values and that all possible hilbert distances
     ;; are seen, plus no two coordinates encode to the same value.

     (define level 4) ; there are (level + 1) actual levels, as 0 is also a level
     (define max-coord (expt 2 (add1 level)))
     (define seen-values (make-hash))
     (for* ([x (in-range max-coord)]
            [y (in-range max-coord)])
       (define d (xy->hilbert-distance level x y))
       (define-values (dx dy) (hilbert-distance->xy level d))
       (check-equal? dx x (format "Failed decoding for x = ~a, y = ~a" x y))
       (check-equal? dy y (format "Failed decoding for x = ~a, y = ~a" x y))
       (check-false (hash-ref seen-values d #f) (format "Duplicate hilbert distance: ~a" d))
       (hash-set! seen-values d #t))

     ;; Check that we have encoded all possible hilbert values.
     (define num-hilbert-distances (expt 4 (add1 level)))
     (for ([d (in-range num-hilbert-distances)])
       (check-true (hash-ref seen-values d #f)
                   (format "Missing hilbert distance: ~a" d))))))

(define pack-unpack-testsuite
  (test-suite
   "Pack/Unpack"
   (test-case "First + Last"
     (check-equal? (first-valid-geoid) 1)
     (check-equal? (last-valid-geoid) #xbfffffffffffffff)
     (check-equal? (sentinel-geoid) (add1 (last-valid-geoid))))
   (test-case "Pack Faces"
     (for ([face (in-range 0 6)])
       (define expected-id (+ (arithmetic-shift face (add1 (* 2 max-level)))
                              (arithmetic-shift 1 (* 2 max-level))))
       (define actual-id (pack face 0 0 max-level))
       (check-equal? expected-id actual-id "Face encoded incorrectly")
       (check-equal? (geoid-face actual-id) face "Face extracted incorrectly")
       (check-equal? (geoid-level actual-id) max-level "Level decoded incorrectly")
       (define-values (face1 x1 y1 level1) (unpack actual-id))
       (check-equal? face1 face "Face decoded incorrectly (1)")
       (check-equal? x1 0 "Bad X coordinate")
       (check-equal? y1 0 "Bad Y coordinate")
       (check-equal? level1 max-level "Bad level")))
   (test-case "Known Encodings"
     (for* ([f (in-range 5)]
            [l (in-range max-level)])
       ;; We know how to calculate the encoding for each face at coordinates
       ;; 0,0, so we can check that pack/unpack produce this value.
       (define expected-encoding
         (+ (arithmetic-shift f (add1 (* 2 max-level)))
            (arithmetic-shift 1 (* 2 l))))
       
       (define p (pack f 0 0 l))
       (check-equal? p expected-encoding "bad encoding")

       (define-values (face x y level) (unpack p))
       (check-equal? face f "bad face")
       (check-equal? x 0 "bad x coordinate")
       (check-equal? y 0 "bad y coordinate")
       (check-equal? level l "bad level")

       ;; Also check geoid-level and geoid-face
       (check-equal? (geoid-level p) l "geoid-level: bad level")
       (check-equal? (geoid-face p) f "geoid-face: bad face")
       (check-equal? (geoid-stride p) (arithmetic-shift 1 (add1 (* 2 l))) "geoid-stride: bad stride")

       (define-values (face1 x1 y1 level1) (unpack (+ p (geoid-stride p))))
       (check-equal? face1 f "bad face")
       (check-equal? level1 l "bad level")
       
       ))
   (test-case "Encodings"
     ;; Check the encodings of the corners for each face at each level pack
     ;; and unpack to the same values
     (for* ([f (in-range 5)]
            [l (in-range max-level)]
            [x (in-list '(0 1))]
            [y (in-list '(0 1))])
       (define max-coord (sub1 (arithmetic-shift 1 (- max-level l))))
       (define id (pack f (* x max-coord) (* y max-coord) l))
       (define-values (face x^ y^ level) (unpack id))
       (check-equal? face f "bad face")
       (check-equal? x^ (* x max-coord) "bad x coordinate")
       (check-equal? y^ (* y max-coord) "bad y coordinate")
       (check-equal? l level "bad level")

       ;; Also check geoid-level and geoid-face here
       (check-equal? (geoid-level id) l "geoid-level: bad level")
       (check-equal? (geoid-face id) f "geoid-face: bad face")
       (check-equal? (geoid-stride id) (arithmetic-shift 1 (add1 (* 2 l))) "geoid-stride: bad stride")

       (let ([next-id (+ id (geoid-stride id))])
         (if (and (= f 5) (= x 1) (= y 1))
             ;; We are on the last ID at this level -- next id should be invalid
             (check-false (valid-geoid? next-id) "next id should not be valid at the end")
             (let-values ([(face1 x1 y1 level1) (unpack next-id)])
               (if (and (= x 1) (= y 0))
                   (check-equal? face1 (add1 f) "geoid-stride: bad face (add1)")
                   (check-equal? face1 f "geoid-stride: bad face"))
               (check-equal? level1 l "geoid-stride: bad level"))))
              
       ))))


(define faces-test-suite
  (test-suite
   "Faces"
   (test-case "Unit Plane Distances"

     ;; The "North Pole"
     (check-= 1.0 ((face-distance top-face) 0 0 1) 1e-10)
     (check-= -1.0 ((face-distance bottom-face) 0 0 1) 1e-10)
     (check-equal? +inf.0 ((face-distance left-face) 0 0 1))
     (check-equal? +inf.0 ((face-distance right-face) 0 0 1))
     (check-equal? +inf.0 ((face-distance front-face) 0 0 1))
     (check-equal? +inf.0 ((face-distance back-face) 0 0 1))

     ;; The "South Pole"
     (check-= -1.0 ((face-distance top-face) 0 0 -1) 1e-10)
     (check-= 1.0 ((face-distance bottom-face) 0 0 -1) 1e-10)
     (check-equal? +inf.0 ((face-distance left-face) 0 0 -1))
     (check-equal? +inf.0 ((face-distance right-face) 0 0 -1))
     (check-equal? +inf.0 ((face-distance front-face) 0 0 -1))
     (check-equal? +inf.0 ((face-distance back-face) 0 0 -1))

     (check-= -1.0 ((face-distance left-face) 0 1 0) 1e-10)
     (check-= 1.0 ((face-distance right-face) 0 1 0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) 0 1 0))
     (check-equal? +inf.0 ((face-distance bottom-face) 0 1 0))
     (check-equal? +inf.0 ((face-distance front-face) 0 1 0))
     (check-equal? +inf.0 ((face-distance back-face) 0 1 0))

     (check-= 1.0 ((face-distance left-face) 0 -1 0) 1e-10)
     (check-= -1.0 ((face-distance right-face) 0 -1 0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) 0 -1 0))
     (check-equal? +inf.0 ((face-distance bottom-face) 0 -1 0))
     (check-equal? +inf.0 ((face-distance front-face) 0 -1 0))
     (check-equal? +inf.0 ((face-distance back-face) 0 -1 0))

     (check-= 1.0 ((face-distance back-face) -1 0 0) 1e-10)
     (check-= -1.0 ((face-distance front-face) -1 0 0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) -1 0 0))
     (check-equal? +inf.0 ((face-distance bottom-face) -1 0 0))
     (check-equal? +inf.0 ((face-distance left-face) -1 0 0))
     (check-equal? +inf.0 ((face-distance right-face) -1 0 0))

     (check-= -1.0 ((face-distance back-face) 1 0 0) 1e-10)
     (check-= 1.0 ((face-distance front-face) 1 0 0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) 1 0 0))
     (check-equal? +inf.0 ((face-distance bottom-face) 1 0 0))
     (check-equal? +inf.0 ((face-distance left-face) 1 0 0))
     (check-equal? +inf.0 ((face-distance right-face) 1 0 0))

     )

   (test-case "Face Ordering"
     ;; Check that the faces in all-faces are in the order of their id, the
     ;; code relies on that (but does not verify it)
     (for ([f (in-list all-faces)]
           [expected-id (in-naturals)])
       (check-equal? (face-id f) expected-id)))

   (test-case "Find Face / Basic"
     (let-values ([(distance face) (find-face 0 0 1)])
       (check-equal? (face-id top-face) (face-id face)))
     (let-values ([(distance face) (find-face 0 0 -1)])
       (check-equal? (face-id bottom-face) (face-id face)))
     (let-values ([(distance face) (find-face 1 0 0)])
       (check-equal? (face-id front-face) (face-id face)))
     (let-values ([(distance face) (find-face -1 0 0)])
       (check-equal? (face-id back-face) (face-id face)))
     (let-values ([(distance face) (find-face 0 1 0)])
       (check-equal? (face-id right-face) (face-id face)))
     (let-values ([(distance face) (find-face 0 -1 0)])
       (check-equal? (face-id left-face) (face-id face))))

   (test-case "Face Project/Unproject"
     ;; Check that the project - unproject functions for each face produce the
     ;; same values. We start with projected values, unproject them, than
     ;; project them again, expecting to get back where we stated.
     ;;
     ;; Note that the unproject are tested by the "Face Continuity" test
     ;; above, so we have greater faith in them.  If this test fails, but the
     ;; "Face Continuity" test does not, the problem is likely with the
     ;; project function.
     (for* ([id (in-range 5)]
            [u (in-range 0 1 0.1)]
            [v (in-range 0 1 0.1)])
       (define face (list-ref all-faces id))
       (define-values (x y z) ((face-uv->xyz face) u v))
       (define d ((face-distance face) x y z))
       (define-values (nu nv) ((face-xyzd->uv face) x y z d))
       (check-= u nu 1e-10 (format "Bad projection for face ~a (x value)" id))
       (check-= v nv 1e-10 (format "Bad projection for face ~a (y value)" id))))

   (test-case "Lat/Lon Encoding"
     ;; This is a degenerate case, all longitudes at the poles will encode to
     ;; the same number
     (define north-pole-geoid (lat-lng->geoid 90 0))
     (for ([lon (in-range -180.0 180.0 2.0)])
       (check-equal? north-pole-geoid (lat-lng->geoid 90 lon)))
     (define south-pole-geoid (lat-lng->geoid -90 0))
     (for ([lon (in-range -180.0 180.0 2.0)])
       (check-equal? south-pole-geoid (lat-lng->geoid -90 lon)))

     ;; The -180.0 and 180.0 are ambiguities and are both encoded as the same
     ;; value, when decoding, we prefer 180.0
     (for ([lat (in-range -89.9 89.9 2.0)])
       (define id1 (lat-lng->geoid lat 180.0))
       (define id2 (lat-lng->geoid lat -180.0))
       (check-equal? id1 id2)
       (define-values (nlat nlon) (geoid->lat-lng id1))
       (check-= nlat lat 1e-5)
       (check-= nlon 180.0 1e-5 (format "failed for lat: ~a" lat)))

     (define latitude-error empty-statistics)
     (define longitude-error empty-statistics)
     (define distance-error empty-statistics)

     ;; Rest of the encodings
     (for* ([lat (in-range -89.9 89.9 2.0)]
            [lon (in-range -179.9 179.9 2.0)])
       (define id (lat-lng->geoid lat lon))
       (define-values (nlat nlon) (geoid->lat-lng id))
       (define dist (map-distance/degrees lat lon nlat nlon))
       #;(printf "lat/lon ~a ~a => ~a => ~a ~a ~%" lat lon id nlat nlon)
       (check < dist 8e-3) ; 6.9 mm is the max error we see (see statistics below)
       (set! latitude-error (update-statistics latitude-error (abs (- lat nlat))))
       (set! longitude-error (update-statistics longitude-error (abs (- lon nlon))))
       (set! distance-error (update-statistics distance-error dist)))

     (pp-statistics latitude-error "Latitude Error")
     (pp-statistics longitude-error "Longitude Error")
     (pp-statistics distance-error "Distance Error (meters)"))

   (test-case "Face Continuity"
     ;; Check that the last hilbert value on one face is next to the first
     ;; hilbert distance on the next face, as well as the GPS points
     ;; corresponding to these are actually neighbours -- this ensures that
     ;; all IDs are continuous across the faces.
     (define max-x (sub1 (expt 2 max-level)))
     ;; NOTE that the last point on the hilbert curve is at (MAX-X, 0), as the
     ;; curve is U shaped.
     (for ([id (in-range 4)]) ; note: up to 4 only, 5 is not connected back to 0.
       (define last-geoid (pack id max-x 0 0))
       (define first-geoid (pack (add1 id) 0 0 0))
       (check-equal?
        (- first-geoid last-geoid) 2
        (format "Discontinuity: last id on face ~a is ~a, first id on face ~a is ~a"
                id last-geoid (add1 id) first-geoid))
       ;; Check that the points are also nearby
       (define-values (last-lat last-lon) (decode id max-x 0 0))
       (define-values (first-lat first-lon) (decode (add1 id) 0 0 0))
       (check-= last-lat first-lat 1e-5
                (format "Latitude discontinuity at face ids ~a / ~a" id (add1 id)))
       (check-= last-lon first-lon 1e-5
                (format "Longitude discontinuity at face ids ~a / ~a" id (add1 id)))))))

(define geoid-manipulation-test-suite
  (test-suite
   "Geoid Manipulation"
   (test-case "leaf-geoid?"
     (for ([k (in-range 10000)])
       (define id (random-geoid 0))
       (check-true (leaf-geoid? id)))
     (for ([k (in-range 10000)])
       (define id (random-geoid (add1 (random 29))))
       (check-false (leaf-geoid? id))))

   (test-case "enclosing-geoid + split-geoid"

     (check-exn exn:fail?
                (lambda () (split-geoid (random-geoid 0))))

     (check-exn exn:fail?
                (lambda () (enclosing-geoid (random-geoid 30))))
     
     (for ([k (in-range 10000)])
       (define id (random-geoid (add1 (random 29))))
       ;; NOTE that split-geoid does not return the four splits in geoid
       ;; order, so we need to sort it if we wish to check the stride.
       (match-define (list g0 g1 g2 g3) (sort (split-geoid id) <))
       (define stride (geoid-stride g0))
       (define level (geoid-level id))
       (check-equal? (geoid-level g0) (sub1 level))
       (check-equal? (geoid-level g1) (sub1 level))
       (check-equal? (geoid-level g2) (sub1 level))
       (check-equal? (geoid-level g3) (sub1 level))
       (check-equal? g1 (+ g0 stride))
       (check-equal? g2 (+ g1 stride))
       (check-equal? g3 (+ g2 stride))
       (check-equal? id (enclosing-geoid g0 (add1 (geoid-level g0))))
       (check-equal? id (enclosing-geoid g1 (add1 (geoid-level g1))))
       (check-equal? id (enclosing-geoid g2 (add1 (geoid-level g2))))
       (check-equal? id (enclosing-geoid g3 (add1 (geoid-level g3))))))

   (test-case "leaf-span"
     (for ([k (in-range 10000)])
       (define id (random-geoid (add1 (random 29))))
       (define level (geoid-level id))
       (define-values (start end) (leaf-span id))
       ;; (printf "~%>>>> id = ~a, start = ~a, end = ~a, level = ~a~%" id start end level)
       (check-true (leaf-geoid? start))
       (check-true
        (or (equal? end (sentinel-geoid))
            (leaf-geoid? (- end (geoid-stride start)))))
       (check-pred < start end)
       (define count
         (if (equal? end (sentinel-geoid))
             (/ (+ (- (sub1 end) start) (geoid-stride start)) (geoid-stride start))
             (/ (- end start) (geoid-stride start))))
       (check-equal? count (expt 2 (* 2 level)))
       (check-equal? (enclosing-geoid start level) id)
       (define last-geoid (if (equal? end (sentinel-geoid))
                              (sub1 end)
                              (- end (geoid-stride start))))
       (check-equal? (enclosing-geoid last-geoid level) id)
       (when (valid-geoid? end)
         (check-not-equal? (enclosing-geoid end level) id))
       (define before (- start (geoid-stride start)))
       (when (valid-geoid? before)
         (check-false (equal? (enclosing-geoid before level) id)))

       ))

   ))

(module+ test
  (require al2-test-runner)
  (parameterize ([current-pseudo-random-generator (make-pseudo-random-generator)])
    (random-seed 42)                    ; ensure our test are deterministic
    (run-tests #:package "geoid"
               #:results-file "test-results-geoid.xml"
               ;;#:exclude '(("Faces"))
               ;; #:only '(("Geoid Manipulation" "leaf-span"))
               hilbert-distance-test-suite
               pack-unpack-testsuite
               faces-test-suite
               geoid-manipulation-test-suite)))
