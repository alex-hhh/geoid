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

(require math/distributions
         math/flonum
         math/statistics
         racket/runtime-path
         rackunit
         al2-test-runner
         (except-in "geodesy.rkt"
                    earth-radius)
         "geoid.rkt"
         "tiling.rkt"
         "vmath.rkt"
         "waypoint-alignment.rkt")

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
     (check-= 1.0 ((face-distance top-face) 0.0 0.0 1.0) 1e-10)
     (check-= -1.0 ((face-distance bottom-face) 0.0 0.0 1.0) 1e-10)
     (check-equal? +inf.0 ((face-distance left-face) 0.0 0.0 1.0))
     (check-equal? +inf.0 ((face-distance right-face) 0.0 0.0 1.0))
     (check-equal? +inf.0 ((face-distance front-face) 0.0 0.0 1.0))
     (check-equal? +inf.0 ((face-distance back-face) 0.0 0.0 1.0))

     ;; The "South Pole"
     (check-= -1.0 ((face-distance top-face) 0.0 0.0 -1.0) 1e-10)
     (check-= 1.0 ((face-distance bottom-face) 0.0 0.0 -1.0) 1e-10)
     (check-equal? +inf.0 ((face-distance left-face) 0.0 0.0 -1.0))
     (check-equal? +inf.0 ((face-distance right-face) 0.0 0.0 -1.0))
     (check-equal? +inf.0 ((face-distance front-face) 0.0 0.0 -1.0))
     (check-equal? +inf.0 ((face-distance back-face) 0.0 0.0 -1.0))

     (check-= -1.0 ((face-distance left-face) 0.0 1.0 0.0) 1e-10)
     (check-= 1.0 ((face-distance right-face) 0.0 1.0 0.0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) 0.0 1.0 0.0))
     (check-equal? +inf.0 ((face-distance bottom-face) 0.0 1.0 0.0))
     (check-equal? +inf.0 ((face-distance front-face) 0.0 1.0 0.0))
     (check-equal? +inf.0 ((face-distance back-face) 0.0 1.0 0.0))

     (check-= 1.0 ((face-distance left-face) 0.0 -1.0 0.0) 1e-10)
     (check-= -1.0 ((face-distance right-face) 0.0 -1.0 0.0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) 0.0 -1.0 0.0))
     (check-equal? +inf.0 ((face-distance bottom-face) 0.0 -1.0 0.0))
     (check-equal? +inf.0 ((face-distance front-face) 0.0 -1.0 0.0))
     (check-equal? +inf.0 ((face-distance back-face) 0.0 -1.0 0.0))

     (check-= 1.0 ((face-distance back-face) -1.0 0.0 0.0) 1e-10)
     (check-= -1.0 ((face-distance front-face) -1.0 0.0 0.0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) -1.0 0.0 0.0))
     (check-equal? +inf.0 ((face-distance bottom-face) -1.0 0.0 0.0))
     (check-equal? +inf.0 ((face-distance left-face) -1.0 0.0 0.0))
     (check-equal? +inf.0 ((face-distance right-face) -1.0 0.0 0.0))

     (check-= -1.0 ((face-distance back-face) 1.0 0.0 0.0) 1e-10)
     (check-= 1.0 ((face-distance front-face) 1.0 0.0 0.0) 1e-10)
     (check-equal? +inf.0 ((face-distance top-face) 1.0 0.0 0.0))
     (check-equal? +inf.0 ((face-distance bottom-face) 1.0 0.0 0.0))
     (check-equal? +inf.0 ((face-distance left-face) 1.0 0.0 0.0))
     (check-equal? +inf.0 ((face-distance right-face) 1.0 0.0 0.0))

     )

   (test-case "Face Ordering"
     ;; Check that the faces in all-faces are in the order of their id, the
     ;; code relies on that (but does not verify it)
     (for ([f (in-list all-faces)]
           [expected-id (in-naturals)])
       (check-equal? (face-id f) expected-id)))

   (test-case "Find Face / Basic"
     (let-values ([(distance face) (find-face 0.0 0.0 1.0)])
       (check-equal? (face-id top-face) (face-id face)))
     (let-values ([(distance face) (find-face 0.0 0.0 -1.0)])
       (check-equal? (face-id bottom-face) (face-id face)))
     (let-values ([(distance face) (find-face 1.0 0.0 0.0)])
       (check-equal? (face-id front-face) (face-id face)))
     (let-values ([(distance face) (find-face -1.0 0.0 0.0)])
       (check-equal? (face-id back-face) (face-id face)))
     (let-values ([(distance face) (find-face 0.0 1.0 0.0)])
       (check-equal? (face-id right-face) (face-id face)))
     (let-values ([(distance face) (find-face 0.0 -1.0 0.0)])
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
            [u (in-range 0.0 1.0 0.1)]
            [v (in-range 0.0 1.0 0.1)])
       (define face (list-ref all-faces id))
       (define-values (x y z) ((face-uv->xyz face) u v))
       (define d ((face-distance face) x y z))
       (define-values (nu nv) ((face-xyzd->uv face) x y z d))
       (check-= u nu 1e-10 (format "Bad projection for face ~a (x value)" id))
       (check-= v nv 1e-10 (format "Bad projection for face ~a (y value)" id))))

   (test-case "Lat/Lon Encoding"
     ;; NOTE that the north pole, south pole and 180 meridians map in the
     ;; middle of the face grid, "between" two adjacent geoids (this is true
     ;; at every split level).  As such, these locations don't have an exact
     ;; geoid and small errors in the input parameters will make it land
     ;; either to the left or to the right of the line.  Also note that `(sin
     ;; pi)` is not exactly zero and `(cos (/ pi 2))` is not exactly zero
     ;; either and this introduces errors.

     ;; This is a degenerate case, all longitudes at the poles will encode to
     ;; the same number
     (define north-pole-geoid (lat-lng->geoid 90.0 0.0))
     (for ([lon (in-range -180.0 180.0 2.0)])
       (check-equal? north-pole-geoid (lat-lng->geoid 90.0 lon)))
     #;(define south-pole-geoid (lat-lng->geoid -90.0 0.0))
     #;(for ([lon (in-range -180.0 180.0 2.0)])
       (check-equal? south-pole-geoid (lat-lng->geoid -90 lon)))

     ;; The -180.0 and 180.0 are ambiguities and are both encoded as the same
     ;; value, when decoding, we prefer 180.0
     (for ([lat (in-range -89.9 89.9 2.0)])
       (define id1 (lat-lng->geoid lat 180.0))
       (define id2 (lat-lng->geoid lat -180.0))

       ;; NOTE: the two geoids are not always equal since the 180 meridian
       ;; falls at the edge of the grid squares and small errors will place
       ;; either to the left or to the right grid cell.  These errors are
       ;; important for tiling, so we keep it this way.

       #;(check-equal? id1 id2)

       (define-values (nlat nlon) (geoid->lat-lng id1))
       (check-= nlat lat 1e-5)
       (check-= (abs nlon) 180.0 5e-5 (format "failed for lat: ~a" lat))

       (define-values (mlat mlon) (geoid->lat-lng id2))
       (check-= mlat lat 1e-5)
       (check-= (abs mlon) 180.0 5e-5 (format "failed for lat: ~a" lat)))

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
       (check < dist 1e-2) ; 9.75 mm is the max error we see (see statistics below)
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
       (define-values (last-lat last-lon)
         (let-values ([(x y z) (decode id max-x 0 0)])
           (unit-vector->lat-lng x y z)))
       (define-values (first-lat first-lon)
         (let-values ([(x y z) (decode (add1 id) 0 0 0)])
           (unit-vector->lat-lng x y z)))
       (check-= last-lat first-lat 1e-5
                (format "Latitude discontinuity at face ids ~a / ~a" id (add1 id)))
       (check-= last-lon first-lon 1e-5
                (format "Longitude discontinuity at face ids ~a / ~a" id (add1 id)))))))

(define geoid-manipulation-test-suite
  (test-suite
   "Geoid Manipulation"

   (test-case "random-geoid"
     (check-exn exn:fail?
                (lambda ()
                  (random-geoid 0 #:parent (random-geoid 0))))
     (check-exn exn:fail?
                (lambda ()
                  (random-geoid 4 #:parent (random-geoid 2))))
     (for ([k (in-range 100)])
       (define plevel (+ 1 (random 29)))
       (define p (random-geoid plevel))
       (for ([j (in-range 100)])
         (define clevel (random plevel))
         (define c (random-geoid clevel #:parent p))
         (check-equal? p (enclosing-geoid c plevel))
         (check-true (contains-geoid? p c)))))

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
       (check < start end)
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

   (test-case "leaf-span*"
     (check-pred null? (leaf-span* '())) ; degenerate case

     (define id (random-geoid (add1 (random 29))))
     (define stride (geoid-stride id))
     (define id1 (+ id stride))
     (define id2 (+ id1 stride))

     (define span1 (leaf-span* (list id id2 id1))) ; NOTE: not in order!
     (check-equal? (length span1) 1)
     (match-define (cons s e) (car span1))
     (check-equal? s (let-values ([(b e) (leaf-span id)]) b))
     (check-equal? e (let-values ([(b e) (leaf-span id2)]) e))

     (define span2 (leaf-span* (list id2 id)))
     (check-equal? (length span2) 2)
     (match-define (list (cons s1 e1) (cons s2 e2)) span2)
     (let-values (([b e] (leaf-span id)))
       (check-equal? s1 b)
       (check-equal? e1 e))
     (let-values (([b e] (leaf-span id2)))
       (check-equal? s2 b)
       (check-equal? e2 e)))

   ))

(define (median x1 y1 z1 x2 y2 z2)
  ;; (printf "median: ~a ~a ~a, ~a ~a ~a~%" x1 y1 z1 x2 y2 z2)
  (define x (+ x1 x2))
  (define y (+ y1 y2))
  (define z (+ z1 z2))
  (define l (sqrt (+ (* x x) (* y y) (* z z))))
  (values (/ x l) (/ y l) (/ z l)))

(define (verify-adjacency geoid neighbours)
  (define-values (x y z) (geoid->unit-vector geoid))
  (define level (geoid-level geoid))
  (for ([n (in-list neighbours)])
    (define-values (nx ny nz) (geoid->unit-vector n))
    (define-values (mx my mz) (median x y z nx ny nz))
    (define-values (face ix iy _level) (encode mx my mz level))
    (define g (pack face ix iy level))
    (check-true
     (or (equal? g geoid)
         (not (equal? #f (member g neighbours))))
     (format "bad neighbour for ~a: ~a" geoid n))))

(define adjacent-test-suite
  (test-suite
   "adjacent-geoids"
   (test-case "faces"
     (for ([face (in-range 6)])
       (define geoid (pack face 0 0 max-level))
       (define opposite (opposite-face face))
       (define neighbours (adjacent-geoids geoid))
       ;; Face level geoids have four unique adjacent geoids
       (check-equal? (length (remove-duplicates neighbours)) 4)
       (verify-adjacency geoid neighbours)))

   (test-case "corners"
     ;; NOTE: at level 30 we only have faces.
     (for* ([level (in-range max-level)]
            [face (in-range 6)]
            [corner (in-list '((0 0) (0 1) (1 1) (1 0)))])
       (match-define (list x y) corner)
       (define m (sub1 (level-info-max-coord (vector-ref level-information level))))
       (define geoid (pack face (* m x) (* m y) level))
       (define neighbours (adjacent-geoids geoid))
       ;; Corners have seven unique adjacent geoids
       (check-equal? (length (remove-duplicates neighbours)) 7)
       (verify-adjacency geoid neighbours)))

   (test-case "edges"
     ;; NOTE: at level 29, there are 4 geoids in a face, all of them are
     ;; corners.
     (for* ([level (in-range (sub1 max-level))]
            [face (in-range 6)])

       (define max-coordinate
         (level-info-max-coord (vector-ref level-information level)))
       (define test-coordinate
         (exact-truncate (/ max-coordinate 2)))
       ;; NOTE: all geoids on the edges should have 8 neighbours

       (define top-edge-geoid (pack face test-coordinate (sub1 max-coordinate) level))
       (define top-neighbours (adjacent-geoids top-edge-geoid))
       (check-equal? (length (remove-duplicates top-neighbours)) 8)
       (verify-adjacency top-edge-geoid top-neighbours)

       (define bottom-edge-geoid (pack face test-coordinate 0 level))
       (define bottom-neigbours (adjacent-geoids bottom-edge-geoid))
       (check-equal? (length (remove-duplicates bottom-neigbours)) 8)
       (verify-adjacency bottom-edge-geoid bottom-neigbours)

       (define left-edge-geoid (pack face 0 test-coordinate level))
       (define left-neigbours (adjacent-geoids left-edge-geoid))
       (check-equal? (length (remove-duplicates left-neigbours)) 8)
       (verify-adjacency left-edge-geoid left-neigbours)

       (define right-edge-geoid (pack face (sub1 max-coordinate) test-coordinate level))
       (define right-neigbours (adjacent-geoids right-edge-geoid))
       (check-equal? (length (remove-duplicates right-neigbours)) 8)
       (verify-adjacency right-edge-geoid right-neigbours)

       ))

   (test-case "middle"
     (for* ([level (in-range (sub1 max-level))]
            [face (in-range 6)])

       (define max-coordinate
         (level-info-max-coord (vector-ref level-information level)))
       (define test-coordinate
         (exact-truncate (/ max-coordinate 2)))
       ;; NOTE: all geoids in the middle of the face should have 8 neighbours

       (define geoid (pack face test-coordinate test-coordinate level))
       (define neigbours (adjacent-geoids geoid))
       (check-equal? (length (remove-duplicates neigbours)) 8)
       (verify-adjacency geoid neigbours)))

   (test-case "random"
     (for* ([level (in-range (add1 max-level))]
            [k (in-range 100)])
       (define geoid (random-geoid level))
       (define neigbours (adjacent-geoids geoid))
       (verify-adjacency geoid neigbours)))

   ))

;; Dams Challenge, Final Climb (Crossing Brookton Highway)
;;
;; s1 -- GPX track from Strava
;; t1 -- Recorded track
;;
(define-runtime-path file-wa-s1 "./test-data/s1.rktd")
(define-runtime-path file-wa-t1 "./test-data/t1.rktd")

;; Mundaring Weir Climb-out
;;
;; s2 -- GPX track from Strava (Kalamunda 100)
;; t2 -- Recorded track (Kalamunda 100, 2021)
;;
(define-runtime-path file-wa-s2 "./test-data/s2.rktd")
(define-runtime-path file-wa-t2 "./test-data/t2.rktd")

;; Lake Gwelup Small Hill
;;
;;
;; s3 -- Extracted GPX route (first lap)
;; t3 -- Recorded track (second lap)
;;
(define-runtime-path file-wa-s3 "./test-data/s3.rktd")
(define-runtime-path file-wa-t3 "./test-data/t3.rktd")

;; Herdsman Loop partial east
;;
;;
;; s4 -- Extracted GPX route (first lap)
;; t4-a -- Recorded track, correct
;; t4-b -- Recorded track, goes on the west side
;;
(define-runtime-path file-wa-s4 "./test-data/s4.rktd")
(define-runtime-path file-wa-t4a "./test-data/t4-a.rktd")
(define-runtime-path file-wa-t4b "./test-data/t4-b.rktd")

(define data-s1 (call-with-input-file file-wa-s1 read))
(define data-t1 (call-with-input-file file-wa-t1 read))
(define data-s2 (call-with-input-file file-wa-s2 read))
(define data-t2 (call-with-input-file file-wa-t2 read))
(define data-s3 (call-with-input-file file-wa-s3 read))
(define data-t3 (call-with-input-file file-wa-t3 read))
(define data-s4 (call-with-input-file file-wa-s4 read))
(define data-t4a (call-with-input-file file-wa-t4a read))
(define data-t4b (call-with-input-file file-wa-t4b read))

(define waypoint-alignment-test-suite
  (test-suite
   "waypoint-alignment"
   (test-case "dtw self-cost"
     (for ([data (in-list (list data-s1 data-t1 data-s2 data-t2
                                data-s3 data-t3 data-s4 data-t4a data-t4b))])
       (define cost (dtw data data))
       ;; NOTE: the cost of a path alignment against itself should be 0, but
       ;; we have some errors in calculating distances, this error seems to be
       ;; just under 4cm for each waypoint
       (check-= (/ cost (vector-length data)) 0 0.04)))
   (test-case "dtw/memory-efficient self-cost"
     (for ([data (in-list (list data-s1 data-t1 data-s2 data-t2
                                data-s3 data-t3 data-s4 data-t4a data-t4b))])
       (define cost (dtw data data))
       ;; NOTE: the cost of a path alignment against itself should be 0, but
       ;; we have some errors in calculating distances, this error seems to be
       ;; just under 4cm for each waypoint
       (check-= (/ cost (vector-length data)) 0 0.04)))

   (test-case "dtw + dtw/memory-efficient equivalence"
     (define (check-cost s t name)
       (define cost-1 (dtw s t))
       (define cost-2 (dtw/memory-efficient s t))
       ;; NOTE: they should be identical!
       (check-= cost-1 cost-2 1e-5 name))

     (check-cost data-s1 data-t1 "s1/t1")
     (check-cost data-s2 data-t2 "s2/t2")
     (check-cost data-s3 data-t3 "s3/t3")
     (check-cost data-s4 data-t4a "s4/t4a")
     (check-cost data-s4 data-t4b "s4/t4b"))

   ))

(define (pp p c)
  (define-values (face ix iy level) (unpack (cell-id c)))
  (printf "*** ~x (~a, ~a, ~a, ~a)~%" (cell-id c) face ix iy level)
  (printf "    cell: ~a~%" c)
  (printf "    point: ~a~%" p)
  (printf "    dot-product: ~a~%" (map (lambda (v) (dot-product p v)) (cell-edges/raw c))))

(define cells-test-suite
  (test-suite
   "cells"
   (test-case "point-inside-cell?"

     (for ([face (in-range 6)])
       (define g (pack face 0 0 30))
       (define c (make-cell g))
       (define-values (x y z) (geoid->unit-vector g))
       (define p (flvector x y z))
       (check-true (point-inside-cell? p c)
                   (format "center point should be inside (geoid #x~x)" g))
       (for ([r (in-list (cell-corners c))])
         (check-true (point-inside-cell? r c)
                     (format "corner point should be inside (geoid #x~x)" g)))
       (for ([h (in-list (adjacent-geoids g))])
         (define-values (x y z) (geoid->unit-vector h))
         (check-false (point-inside-cell? (flvector x y z) c)
                      (format "neighbor center should be outside (geoid #x~x, neighbor #x~x)"
                              g h))))

     (for ([_k (in-range 10000)])
       (define level (random 30))
       (define g (random-geoid level))
       (define-values (x y z) (geoid->unit-vector g))
       (define c (make-cell g))
       (check-true (point-inside-cell? (flvector x y z) c)
                   (format "center point should be inside (geoid #x~x)" g))
       ;; This has errors for small geoids...
       (for ([r (in-list (cell-corners c))])
         (check-true (point-inside-cell? r c)
                     (format "corner point should be inside (geoid #x~x) ~a" g r)))
       (for ([h (in-list (adjacent-geoids g))])
         (define-values (x y z) (geoid->unit-vector h))
         (check-true (point-outside-cell? (flvector x y z) c)
                     (format "neighbor center should be outside (geoid #x~x, neighbor #x~x)"
                             g h)))))))

;; A loop track around the swan river in Perth, used to test geoid covering
(define-runtime-path swan-track-file "./test-data/swan-track.rktd")
(define swan-track (call-with-input-file swan-track-file read))
(define-runtime-path swan-track-covering-file "./test-data/swan-track-covering.rktd")
(define swan-track-expected-covering (call-with-input-file swan-track-covering-file read))

;; A small island on the France-Spain border, defined in ClockWise order
;;
;; https://en.wikipedia.org/wiki/Pheasant_Island
(define-runtime-path pheasant-island-file "./test-data/pheasant-island.rktd")
(define pheasant-island (call-with-input-file pheasant-island-file read))

;; A small Italian region inside Switzerland, defined in CounterClockWise
;; order.
;;
;; https://en.wikipedia.org/wiki/Campione_d'Italia
(define-runtime-path campione-ditalia-file "./test-data/campione-ditalia.rktd")
(define campione-ditalia (call-with-input-file campione-ditalia-file read))

;; This is the Antarctica/McMurdo time zone which contains the South Pole.  It
;; is defined such that it looks nice on a Mercator map, but it goes through
;; the South Pole, and uses two separate values for it, (-90, -180) and (-90,
;; 180)
(define-runtime-path antarctica-mcmurdo-file "./test-data/mcmurdo.rktd")
(define antarctica-mcmurdo (call-with-input-file antarctica-mcmurdo-file read))
(define-runtime-path antarctica-mcmurdo-contains-file "./test-data/mcmurdo-contains.rktd")
(define antarctica-mcmurdo-contains (call-with-input-file antarctica-mcmurdo-contains-file read))
(define-runtime-path antarctica-mcmurdo-intersects-file "./test-data/mcmurdo-intersects.rktd")
(define antarctica-mcmurdo-intersects (call-with-input-file antarctica-mcmurdo-intersects-file read))

;; This is the first loop in the Pacific/Funafuti time zone, it has a line
;; along the 180 meridian, which gave me endless trouble with the triangle
;; signs.
(define-runtime-path pacific-funafuti-file "./test-data/funafuti.rktd")
(define pacific-funafuti (call-with-input-file pacific-funafuti-file read))

(define-runtime-path weir-climbout-file "./test-data/weir-climbout.rktd")
(define weir-climbout (call-with-input-file weir-climbout-file read))
(define-runtime-path weir-climbout-covering-file "./test-data/weir-climbout-covering.rktd")
(define weir-climbout-covering (call-with-input-file weir-climbout-covering-file read))

(define (make-circular-track lat lon radius
         #:segments [segments 24] #:closed? [closed? #t])
  (define step (- (/ 360.0 segments)))
  (define track
    (for/list ([bearing (in-range 0.0 -360.0 step)])
      (define-values (clat clon) (destination-point lat lon bearing radius))
      (vector clat clon)))
  (if closed?
      (append track (list (car track)))
      track))

(define tiling-test-suite
  (test-suite
   "tiling"

   (test-case "angle-contains-vertex?"
     (define v0 (flvector -0.3704386151255099 0.7632268034408783 -0.5293959566650591))
     (define v1 (flvector -0.37044282581595456 0.763232968638743 -0.5293841217720663))
     (define v2 (flvector -0.37044282581595456 0.763232968638743 -0.5293841217720663))

     ;; I'm still not sure what it means for an angle to contain its own
     ;; middle vertex, but the property that we are interested in is that, if
     ;; an angle contains a vertex that the reverse angle does not.
     (check-false (angle-contains-vertex? v0 v1 v2))
     (check-true (angle-contains-vertex? v2 v1 v0))


     (define a (flvector 1.0 0.0 0.0))
     (define b (flvector 0.0 1.0 0.0))
     (define ref_b (reference-direction b))
     (check-false (angle-contains-vertex? a b a))
     (check-true (angle-contains-vertex? ref_b b a))
     (check-false (angle-contains-vertex? a b ref_b))

     (define-values (lat lon) (unit-vector->lat-lng 0.0 1.0 0.0))
     (define points (map ->unit-vector (make-circular-track lat lon 100)))
     (define count
       (for/sum ([p1 (in-list points)]
                 [p2 (in-list (cdr points))])
         (define c1 (angle-contains-vertex? p2 b p1))
         (define c2 (angle-contains-vertex? p1 b p2))
         ;; One (and only one) of the angles (P1 -> B -> P2) and (P2 -> B ->
         ;; P1) will contain the vertex B.
         (check xor c1 c2)
         (if c1 1 0)))
     (check-equal? count 1))

   (test-case "segment and vertex crossing"
     (define p0 (->unit-vector #(0.0 5.0)))
     (define p1 (->unit-vector #(0.0 10.0)))
     (define p2 (->unit-vector #(0.0 20.0)))

     (define s0 (make-segment p1 p2))
     (define s1 (make-segment p2 p0))

     (define north-pole (flvector 0.0 0.0 1.0))
     (define south-pole (flvector 0.0 0.0 -1.0))

     ;; The p3 -- north pole segment should cross s0, while the p3 --
     ;; south-pole should not.
     (define p3 (->unit-vector #(-10.0 15.0)))
     (check-equal? (segment-crossing-sign (make-segment p3 north-pole) s0) 1)
     (check-equal? (segment-crossing-sign (make-segment north-pole p3) s0) 1)
     (check-equal? (segment-crossing-sign (make-segment p3 south-pole) s0) -1)
     (check-equal? (segment-crossing-sign (make-segment south-pole p3) s0) -1)

     ;; The line from p4 to the north-pole should touch the start of S0 and
     ;; the end of S1, it is only intersecting S1, i.e. P1 is on S1 but not
     ;; S0.
     (define p4 (->unit-vector #(-10.0 10.0)))
     (check-equal? (segment-crossing-sign (make-segment p4 north-pole) s0) -1)
     (check-equal? (segment-crossing-sign (make-segment north-pole p4) s0) -1)
     (check-equal? (segment-crossing-sign (make-segment p4 north-pole) s1) 1)
     (check-equal? (segment-crossing-sign (make-segment north-pole p4) s1) 1)


     ;; The line from p5 -- north-pole should touch the end of S0
     (define p5 (->unit-vector #(-10.0 20.0)))
     (check-equal? (segment-crossing-sign (make-segment p5 north-pole) s0) -1)
     (check-equal? (segment-crossing-sign (make-segment north-pole p5) s0) -1)
     (check-equal? (segment-crossing-sign s0 (make-segment p5 north-pole)) -1)
     (check-equal? (segment-crossing-sign s0 (make-segment north-pole p5)) -1)

     ;; p6 is on S0, but only one segment is considered to be crossing (the
     ;; north pole one).
     (define p6 (->unit-vector #(0.0 15.0)))
     (check-equal? (segment-crossing-sign (make-segment p6 north-pole) s0) 1)
     (check-equal? (segment-crossing-sign (make-segment north-pole p6) s0) 1)
     (check-equal? (segment-crossing-sign (make-segment p6 south-pole) s0) -1)
     (check-equal? (segment-crossing-sign (make-segment south-pole p6) s0) -1)
     (check-equal? (segment-crossing-sign s0 (make-segment p6 north-pole)) 1)
     (check-equal? (segment-crossing-sign s0 (make-segment north-pole p6)) 1)
     (check-equal? (segment-crossing-sign s0 (make-segment p6 south-pole)) -1)
     (check-equal? (segment-crossing-sign s0 (make-segment south-pole p6)) -1)


     ;; segment-crossing-sign returns 0 if two points from the segments are
     ;; the same...

     (let ([s1 (make-segment p1 north-pole)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0))
       (check-false (vertex-crossing? s0 s1)))

     (let ([s1 (make-segment north-pole p1)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0))
       (check-false (vertex-crossing? s0 s1)))

     (let ([s1 (make-segment p1 south-pole)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-false (vertex-crossing? s1 s0))
       (check-true (vertex-crossing? s0 s1)))

     (let ([s1 (make-segment south-pole p1)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-false (vertex-crossing? s1 s0))
       (check-true (vertex-crossing? s0 s1)))

     ;; p2 is on the end of S0
     (let ([s1 (make-segment p2 north-pole)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0))
       (check-false (vertex-crossing? s0 s1)))

     (let ([s1 (make-segment north-pole p2)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0))
       (check-false (vertex-crossing? s0 s1)))

     (let ([s1 (make-segment p2 south-pole)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0))
       (check-false (vertex-crossing? s0 s1)))

     (let ([s1 (make-segment south-pole p2)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0))
       (check-false (vertex-crossing? s0 s1)))

     ;; A segment will cross itself (even if it is reversed)
     (check-equal? (segment-crossing-sign s0 s0) 0)
     (check-true (vertex-crossing? s0 s0))

     (let ([s1 (make-segment p2 p1)])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0)))

     ;; Segments are co-linear
     (let ([s1 (make-segment p1 (->unit-vector #(0.0 15.0)))])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-false (vertex-crossing? s1 s0))
       (check-true (vertex-crossing? s0 s1)))
     (let ([s1 (make-segment p1 (->unit-vector #(0.0 25.0)))])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-true (vertex-crossing? s1 s0))
       (check-false (vertex-crossing? s0 s1)))
     (let ([s1 (make-segment p1 (->unit-vector #(0.0 -10.0)))])
       (check-equal? (segment-crossing-sign s1 s0) 0)
       (check-false (vertex-crossing? s1 s0))
       (check-true (vertex-crossing? s0 s1))))

   (test-case "is-point-inside-loop?"
     (define test-point (flvector 0.0 0.0 1.0))

     ;; NOTE: the swan track is defined clockwise, so we expect the north pole
     ;; to be inside the track by default...

     (let* ([track swan-track]
            [segments (track->segments track)])
       (define-values (origin is-inside?) (get-suitable-origin segments))
       (check-true
        (is-point-inside-loop? test-point segments origin is-inside?)))

     (let* ([track (reverse swan-track)]
            [segments (track->segments track)])
       (define-values (origin is-inside?) (get-suitable-origin segments))
       (check-false
        (is-point-inside-loop? test-point segments origin is-inside?)))

     (let* ([track '(#(0.0 -10.0) #(-10.0 0.0) #(0.0 10.0) #(10.0 0.0))]
            [segments (track->segments track #t)]
            [origin (flvector 0.0 0.0 1.0)])

       ;; The test segment exits through a corner of the track
       (let ([test-point (->unit-vector #(0.0 0.0))])
         (check-true
          (is-point-inside-loop? test-point segments origin #f)))

       ;; The test segment touches a point of our track
       (let ([test-point (->unit-vector #(-10.0 -10.0))])
         (check-false
          (is-point-inside-loop? test-point segments origin #f)))))

   (test-case "guess-winding-order"
     (check-equal? (guess-winding-order pheasant-island) 'cw)
     (check-equal? (guess-winding-order campione-ditalia) 'ccw))

   (test-case "empty-cap"
     (check-equal? (geoid-tiling-for-region the-empty-cap 13 18) '()))
   (test-case "singular-cap"
     ;; This is a cap containing only its center, tiling it should produce a
     ;; single geoid at the minimum level.
     (define singular-cap (make-spherical-cap -31.96355 115.84529 0))
     (define covering (geoid-tiling-for-region singular-cap 13 18))
     (check-equal? covering '(2989588974810431488)))
   (test-case "non-singular-cap"
     ;; The test case was created by generating a covering, visually
     ;; inspecting it on the map (see the tools.rkt file) and saving it here.
     ;; If an unexpected result is produced here, the results can be
     ;; visualized on a map to determine what is wrong
     (define test-cap (make-spherical-cap -31.96355 115.84529 500))
     (define covering (geoid-tiling-for-region test-cap 13 18))
     (define expected-result
       '(#x297cd8ab44000000 #x297cd8ab30000000 #x297cd8ab10000000 #x297cd8ab6c000000
         #x297cd8ab74000000 #x297cd8aac0000000 #x297cd8aa40000000 #x297cd8ab9c000000
         #x297cd8ab94000000 #x297cd8abb0000000 #x297cd8abd0000000 #x297cd8a97c000000
         #x297cd8a990000000 #x297cd8a9f0000000 #x297d20ac84000000 #x297d20ac9c000000
         #x297d20aca4000000 #x297d20acac000000 #x297d20aa14000000 #x297d20aa0c000000
         #x297d20aa34000000 #x297d20aa3c000000 #x297d20aa50000000 #x297d20aa70000000
         #x297d20abc4000000 #x297d20abcc000000 #x297d20ab90000000 #x297d20abb0000000
         #x297d20ab40000000 #x297d20aac0000000 #x297d275540000000 #x297d2754c0000000
         #x297d275450000000 #x297d275470000000 #x297d275414000000 #x297d275430000000
         #x297d2755c0000000 #x297d275354000000 #x297d27535c000000 #x297d275364000000
         #x297d27537c000000 #x297d275610000000 #x297d275674000000 #x297d27567c000000
         #x297d27566c000000 #x297d275684000000 #x297cdf5590000000 #x297cdf55a4000000
         #x297cdf5540000000 #x297cdf54e4000000 #x297cdf54d0000000))
     (check-equal? (sort covering <) (sort expected-result <)))
   (test-case "open-polyline"
     (define region (make-open-polyline weir-climbout))
     (define covering (geoid-tiling-for-region region 13 18))
     (define expected-result weir-climbout-covering)
     (check-equal? (sort covering <) (sort expected-result <)))
   (test-case "closed-polyline"
     (define expected-result swan-track-expected-covering)

     (define region (make-closed-polyline swan-track #:ccw? #f))
     (define covering (geoid-tiling-for-region region 13 18))
     (check-equal? (sort covering <) (sort expected-result <))

     (define indexed-region (make-closed-polyline
                             swan-track #:ccw? #f #:force-indexed? #t))
     (define indexed-covering (geoid-tiling-for-region indexed-region 13 18))
     (check-equal? (sort indexed-covering <) (sort expected-result <)))

   (test-case "closed-polyline-around-north-pole"
     (define track
       '(#(89.95 0.0) #(89.95 90.0) #(89.95 180.0) #(89.95 -90.0)))
     (define expected-result
       '(#x4555552b00000000 #x4555552d00000000 #x4555552f00000000 #x4555553400000000 #x4555553900000000
         #x4555553b00000000 #x4555553f00000000 #x4555555000000000 #x4555556400000000 #x4555556900000000
         #x4555556b00000000 #x4555556f00000000 #x4555557b00000000 #x4555557d00000000 #x4555557f00000000
         #x4fffff8100000000 #x4fffff8500000000 #x4fffff8700000000 #x4fffff8c00000000 #x4fffff9100000000
         #x4fffff9300000000 #x4fffff9500000000 #x4fffffc100000000 #x4fffffc300000000 #x4fffffc500000000
         #x4fffffd100000000 #x4fffffd500000000 #x4fffffd700000000 #x4fffffdc00000000 #x4ffffff000000000
         #x5000001000000000 #x5000002400000000 #x5000002900000000 #x5000002b00000000 #x5000002f00000000
         #x5000003b00000000 #x5000003d00000000 #x5000003f00000000 #x5000006b00000000 #x5000006d00000000
         #x5000006f00000000 #x5000007400000000 #x5000007900000000 #x5000007b00000000 #x5000007f00000000
         #x5aaaaa8100000000 #x5aaaaa8300000000 #x5aaaaa8500000000 #x5aaaaa9100000000 #x5aaaaa9500000000
         #x5aaaaa9700000000 #x5aaaaa9c00000000 #x5aaaaab000000000 #x5aaaaac100000000 #x5aaaaac500000000
         #x5aaaaac700000000 #x5aaaaacc00000000 #x5aaaaad100000000 #x5aaaaad300000000 #x5aaaaad500000000))
     (define region (make-closed-polyline track))
     (define covering (geoid-tiling-for-region region 16 18))
     (check-equal? (sort covering <) (sort expected-result <))

     (define indexed-region (make-closed-polyline track #:force-indexed? #t))
     (define indexed-covering (geoid-tiling-for-region indexed-region 16 18))
     (check-equal? (sort indexed-covering <) (sort expected-result <)))

   (test-case "Antarctica/McMurdo covering"
     (define region (make-closed-polyline antarctica-mcmurdo #:ccw? #f))
     (define-values (contains intersects) (coarse-geoid-tiling-for-region region 24))
     (define expected-contains antarctica-mcmurdo-contains)
     (define expected-intersects antarctica-mcmurdo-intersects)
     (check-equal? contains expected-contains)
     (check-equal? intersects expected-intersects))

   (test-case "Pacific/Funafuti covering"
     (define region (make-closed-polyline pacific-funafuti))
     (define-values (contains intersects) (coarse-geoid-tiling-for-region region 24))
     (define expected-contains '())
     (define expected-intersects '(#x7027000000000000))
     (check-equal? contains expected-contains)
     (check-equal? intersects expected-intersects))

   (test-case "overlapping union"
     (define cap1 (make-spherical-cap -31.96355 115.84529 500))
     (define cap2 (make-spherical-cap -31.96537 115.84900 500))
     (define u (join-regions cap1 cap2))
     (define covering (geoid-tiling-for-region u 13 18))
     (define expected-result
       '(#x297cd8ab44000000 #x297cd8ab30000000 #x297cd8ab10000000 #x297cd8ab6c000000 #x297cd8ab74000000
         #x297cd8aac0000000 #x297cd8aa40000000 #x297cd8ab9c000000 #x297cd8ab94000000 #x297cd8abb0000000
         #x297cd8abd0000000 #x297cd8a97c000000 #x297cd8a990000000 #x297cd8a9f0000000 #x297d20a914000000
         #x297d20a90c000000 #x297d20a930000000 #x297d20a950000000 #x297d20a970000000 #x297d20a990000000
         #x297d20a9f4000000 #x297d20ae5c000000 #x297d20ae64000000 #x297d20ae74000000 #x297d20ae7c000000
         #x297d20aed0000000 #x297d20aee4000000 #x297d20aef4000000 #x297d20aeec000000 #x297d20ae90000000
         #x297d20aeb0000000 #x297d20ac40000000 #x297d20add0000000 #x297d20adb0000000 #x297d20ad10000000
         #x297d20ad30000000 #x297d20acc0000000 #x297d20ab00000000 #x297d275540000000 #x297d2754c0000000
         #x297d275450000000 #x297d275470000000 #x297d275414000000 #x297d275430000000 #x297d2755c0000000
         #x297d275350000000 #x297d275330000000 #x297d275314000000 #x297d275370000000 #x297d2752d4000000
         #x297d275384000000 #x297d275610000000 #x297d275674000000 #x297d27567c000000 #x297d27566c000000
         #x297d275684000000 #x297cdf5590000000 #x297cdf55a4000000 #x297cdf5540000000 #x297cdf54e4000000
         #x297cdf54d0000000))
     (check-equal? (sort covering <) (sort expected-result <)))

   (test-case "disjoint union"
     (define cap1 (make-spherical-cap -31.96355 115.84529 500))
     (define cap3 (make-spherical-cap -31.97429 115.86697 500))
     (define u (join-regions cap1 cap3))
     (define covering (geoid-tiling-for-region u 13 18))
     (define expected-result
       '(#x297cd8ab44000000 #x297cd8ab30000000 #x297cd8ab10000000 #x297cd8ab6c000000 #x297cd8ab74000000
         #x297cd8aac0000000 #x297cd8aa40000000 #x297cd8ab9c000000 #x297cd8ab94000000 #x297cd8abb0000000
         #x297cd8abd0000000 #x297cd8a97c000000 #x297cd8a990000000 #x297cd8a9f0000000 #x297d20be2c000000
         #x297d20be34000000 #x297d20be3c000000 #x297d20be50000000 #x297d20be6c000000 #x297d20be64000000
         #x297d20bfc0000000 #x297d20bf40000000 #x297d20bec0000000 #x297d20b940000000 #x297d20b8c0000000
         #x297d20b850000000 #x297d20b870000000 #x297d20b814000000 #x297d20b80c000000 #x297d20b834000000
         #x297d20b83c000000 #x297d20b99c000000 #x297d20b9b0000000 #x297d20b9cc000000 #x297d20ac84000000
         #x297d20ac9c000000 #x297d20aca4000000 #x297d20acac000000 #x297d20aa14000000 #x297d20aa0c000000
         #x297d20aa34000000 #x297d20aa3c000000 #x297d20aa50000000 #x297d20aa70000000 #x297d20abc4000000
         #x297d20abcc000000 #x297d20ab90000000 #x297d20abb0000000 #x297d20ab40000000 #x297d20aac0000000
         #x297d20ea0c000000 #x297d20ea70000000 #x297d20eac4000000 #x297d20eaec000000 #x297d20ea90000000
         #x297d20eab0000000 #x297d275540000000 #x297d2754c0000000 #x297d275450000000 #x297d275470000000
         #x297d275414000000 #x297d275430000000 #x297d2755c0000000 #x297d275354000000 #x297d27535c000000
         #x297d275364000000 #x297d27537c000000 #x297d275610000000 #x297d275674000000 #x297d27567c000000
         #x297d27566c000000 #x297d275684000000 #x297cdf5590000000 #x297cdf55a4000000 #x297cdf5540000000
         #x297cdf54e4000000 #x297cdf54d0000000 #x297d20c040000000 #x297d20c1d0000000 #x297d20c1e4000000
         #x297d20c1f4000000 #x297d20c1ec000000 #x297d20c190000000 #x297d20c1b0000000 #x297d20c110000000
         #x297d20c170000000 #x297d20c144000000 #x297d20c15c000000 #x297d20c14c000000 #x297d20c130000000
         #x297d20c0c0000000 #x297d20c740000000 #x297d20c6d0000000 #x297d20c6b4000000 #x297d20c6bc000000
         #x297d20c6e4000000 #x297d20c6ec000000 #x297d20c6fc000000 #x297d20c790000000 #x297d20c7a4000000
         #x297d20c7ac000000 #x297d20c7bc000000 #x297d20c7f4000000 #x297d209550000000 #x297d209564000000
         #x297d20957c000000))
     (check-equal? (sort covering <) (sort expected-result <)))

   (test-case "intersection"
     (define cap1 (make-spherical-cap -31.96355 115.84529 500))
     (define cap2 (make-spherical-cap -31.96537 115.84900 500))
     (define u (intersect-regions cap1 cap2))
     (define covering (geoid-tiling-for-region u 13 18))
     (define expected-result
       '(#x297cd8aaac000000 #x297cd8aaa4000000 #x297cd8aa9c000000 #x297d20ac84000000 #x297d20ac9c000000
         #x297d20aca4000000 #x297d20acac000000 #x297d20aa14000000 #x297d20aa0c000000 #x297d20aa34000000
         #x297d20aa3c000000 #x297d20aa50000000 #x297d20aa70000000 #x297d20abc4000000 #x297d20abcc000000
         #x297d20ab90000000 #x297d20abb0000000 #x297d20ab40000000 #x297d20aac0000000 #x297d275540000000
         #x297d2754c0000000 #x297d275450000000 #x297d275464000000 #x297d27547c000000 #x297d27559c000000
         #x297d2755a4000000 #x297d2755ac000000 #x297d2755b4000000 #x297d275354000000 #x297d27535c000000
         #x297d275364000000 #x297d27537c000000 #x297cdf558c000000 #x297cdf5584000000 #x297cdf557c000000
         #x297cdf5564000000 #x297cdf555c000000 #x297cdf5554000000))
     (check-equal? (sort covering <) (sort expected-result <)))

   (test-case "disjoint intersection"
     (define cap1 (make-spherical-cap -31.96355 115.84529 500))
     (define cap3 (make-spherical-cap -31.97429 115.86697 500))
     (define u (intersect-regions cap1 cap3))
     (define covering (geoid-tiling-for-region u 13 18))
     ;; Caps are disjoint, expecting an empty intersection set
     (define expected-result '())
     (check-equal? (sort covering <) (sort expected-result <)))

   (test-case "subtraction"
     (define cap1 (make-spherical-cap -31.96355 115.84529 500))
     (define cap2 (make-spherical-cap -31.96537 115.84900 500))
     (define u (subtract-regions cap1 cap2))
     (define covering (geoid-tiling-for-region u 13 18))
     (define expected-result
       '(#x297cd8ab44000000 #x297cd8ab30000000 #x297cd8ab10000000 #x297cd8ab6c000000 #x297cd8ab74000000
         #x297cd8aac0000000 #x297cd8aa40000000 #x297cd8ab9c000000 #x297cd8ab94000000 #x297cd8abb0000000
         #x297cd8abd0000000 #x297cd8a97c000000 #x297cd8a990000000 #x297cd8a9f0000000 #x297d20aa0c000000
         #x297d20aa74000000 #x297d275564000000 #x297d275574000000 #x297d27557c000000 #x297d27545c000000
         #x297d275444000000 #x297d27544c000000 #x297d275470000000 #x297d275414000000 #x297d275430000000
         #x297d2755c0000000 #x297d27537c000000 #x297d275610000000 #x297d275674000000 #x297d27567c000000
         #x297d27566c000000 #x297d275684000000 #x297cdf5590000000 #x297cdf55a4000000 #x297cdf5540000000
         #x297cdf54e4000000 #x297cdf54d0000000))
     (check-equal? (sort covering <) (sort expected-result <)))

   (test-case "disjoint subtraction"
     (define cap1 (make-spherical-cap -31.96355 115.84529 500))
     (define cap3 (make-spherical-cap -31.97429 115.86697 500))
     (define u (subtract-regions cap1 cap3))
     (define covering (geoid-tiling-for-region u 13 18))
     (define expected-result
       '(#x297cd8ab44000000 #x297cd8ab30000000 #x297cd8ab10000000 #x297cd8ab6c000000 #x297cd8ab74000000
         #x297cd8aac0000000 #x297cd8aa40000000 #x297cd8ab9c000000 #x297cd8ab94000000 #x297cd8abb0000000
         #x297cd8abd0000000 #x297cd8a97c000000 #x297cd8a990000000 #x297cd8a9f0000000 #x297d20ac84000000
         #x297d20ac9c000000 #x297d20aca4000000 #x297d20acac000000 #x297d20aa14000000 #x297d20aa0c000000
         #x297d20aa34000000 #x297d20aa3c000000 #x297d20aa50000000 #x297d20aa70000000 #x297d20abc4000000
         #x297d20abcc000000 #x297d20ab90000000 #x297d20abb0000000 #x297d20ab40000000 #x297d20aac0000000
         #x297d275540000000 #x297d2754c0000000 #x297d275450000000 #x297d275470000000 #x297d275414000000
         #x297d275430000000 #x297d2755c0000000 #x297d275354000000 #x297d27535c000000 #x297d275364000000
         #x297d27537c000000 #x297d275610000000 #x297d275674000000 #x297d27567c000000 #x297d27566c000000
         #x297d275684000000 #x297cdf5590000000 #x297cdf55a4000000 #x297cdf5540000000 #x297cdf54e4000000
         #x297cdf54d0000000))
     (check-equal? (sort covering <) (sort expected-result <)))

   ))

(define (test-data-points/vincenty lat1 lon1 lat2 lon2)
  (define-values (distance-1 initial-bearing-1 final-bearing-1)
    (vincenty-inverse lon1 lat1 lon2 lat2))
  (define-values (distance-2 initial-bearing-2 final-bearing-2)
    (vincenty-inverse lon2 lat2 lon1 lat1))

  (define-values (dlon2 dlat2 fb2)
    (vincenty-direct lon1 lat1 initial-bearing-1 distance-1))
  (define-values (dlon1 dlat1 fb1)
    (vincenty-direct lon2 lat2 initial-bearing-2 distance-2))

  (check-= distance-1 distance-2 1e-6)
  (check-= initial-bearing-1 (wrap2π (+ final-bearing-2 pi)) 1e-6)
  (check-= final-bearing-1 (wrap2π (+ initial-bearing-2 pi)) 1e-6)
  (check-= lat1 dlat1 1e-6)
  (check-= lon1 dlon1 1e-6)
  (check-= lat2 dlat2 1e-6)
  (check-= lon2 dlon2 1e-6))

(define (test-distance-and-bearing/spherical lat1 lon1 lat2 lon2)
  (define distance-1 (spherical-distance lon1 lat1 lon2 lat2))
  (define distance-2 (spherical-distance lon2 lat2 lon1 lat1))

  (check-= distance-1 distance-2 1e-6)

  (define initial-bearing-1 (spherical-bearing lon1 lat1 lon2 lat2))
  (define-values (lon3 lat3)
    (spherical-destination lon1 lat1 initial-bearing-1 distance-1))
  (check-= lat2 lat3 1e-6)
  (check-= lon2 lon3 1e-6)

  (define initial-bearing-2 (spherical-bearing lon2 lat2 lon1 lat1))
  (define-values (lon4 lat4)
    (spherical-destination lon2 lat2 initial-bearing-2 distance-2))
  (check-= lat1 lat4 1e-6)
  (check-= lon1 lon4 1e-6))

(define (test-midway-point/spherical lat1 lon1 lat2 lon2)
  (define-values (lonm1 latm1) (spherical-midway-point lon1 lat1 lon2 lat2))
  (define-values (lonm2 latm2) (spherical-midway-point lon2 lat2 lon1 lat1))

  (check-= lonm1 lonm2 1e-6)
  (check-= latm1 latm2 1e-6)

  (define first-half (spherical-distance lon1 lat1 lonm1 latm1))
  (define second-half (spherical-distance lonm1 latm1 lon2 lat2))
  (define full (spherical-distance lon1 lat1 lon2 lat2))

  (check-= first-half second-half 1e-6)
  (check-= (+ first-half second-half) full 1e-6))

(define (test-intermediate-point/spherical lat1 lon1 lat2 lon2)
  (define-values (lon-start lat-start)
    (spherical-intermediate-point lon1 lat1 lon2 lat2 0))

  (check-= lon-start lon1 1e-6)
  (check-= lat-start lat1 1e-6)

  (define-values (lon-end lat-end)
    (spherical-intermediate-point lon1 lat1 lon2 lat2 1))

  (check-= lon-end lon2 1e-6)
  (check-= lat-end lat2 1e-6)

  (define-values (lon-mid lat-mid)
    (spherical-intermediate-point lon1 lat1 lon2 lat2 1/2))
  (define-values (lonm latm) (spherical-midway-point lon1 lat1 lon2 lat2))

  (check-= lon-mid lonm 1e-6)
  (check-= lat-mid latm 1e-6))


(define geodesy-test-suite
  (test-suite
    "Geodesy test suite"
    (test-case "vincenty direct - inverse consistency"

      (define iterations 5000)
      (define lat-distribution (uniform-dist (- (/ pi 2)) (/ pi 2)))
      (define lon-distribution (uniform-dist (- pi) pi))

      (for ([lat1 (sample lat-distribution iterations)]
            [lon1 (sample lon-distribution iterations)]
            [lat2 (sample lat-distribution iterations)]
            [lon2 (sample lon-distribution iterations)])
        (test-data-points/vincenty lat1 lon1 lat2 lon2)))

    (test-case "spherical direct - inverse consistency"
      (define iterations 5000)
      (define lat-distribution (uniform-dist (- (/ pi 2)) (/ pi 2)))
      (define lon-distribution (uniform-dist (- pi) pi))

      (for ([lat1 (sample lat-distribution iterations)]
            [lon1 (sample lon-distribution iterations)]
            [lat2 (sample lat-distribution iterations)]
            [lon2 (sample lon-distribution iterations)])
        (test-distance-and-bearing/spherical lat1 lon1 lat2 lon2)
        (test-midway-point/spherical lat1 lon1 lat2 lon2)
        (test-intermediate-point/spherical lat1 lon1 lat2 lon2)))

  ))

(module+ test
  (require al2-test-runner)
  (parameterize ([current-pseudo-random-generator (make-pseudo-random-generator)])
    (random-seed 42)                    ; ensure our test are deterministic
    (run-tests #:package "geoid"
               #:results-file "test-results-geoid.xml"
               ;;#:exclude '(("tiling" "Antarctica/McMurdo covering"))
               ;;#:only '(("tiling" "angle-contains-vertex?"))
               hilbert-distance-test-suite
               pack-unpack-testsuite
               faces-test-suite
               geoid-manipulation-test-suite
               adjacent-test-suite
               waypoint-alignment-test-suite
               cells-test-suite
               tiling-test-suite
               geodesy-test-suite)))
