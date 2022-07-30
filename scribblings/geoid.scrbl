#lang scribble/manual
@require[@for-label[geoid
                    racket/class
                    geoid/waypoint-alignment
                    geoid/geodesy
                    geoid/tiling
                    racket/contract/base
                    racket/base]]

@title{geoid -- Efficient storage and queriying of geographic data}
@author{Alex Harsányi}

@defmodule[geoid]

The @racket[geoid] library allows representing latidude/longitude coordinates
using 64-bit integers.  The integer values are ordered to preseve locality,
which means that these IDs can be stored in databases and used to query and
efficiently retrieve locations which are in a geographic region.

This library is inspired by the @hyperlink["http://s2geometry.io/"]{S2
Geometry} library, in the sense that the integer IDs are determined in the
same way.  Until this libary has its own introduction section, you can read
@hyperlink["http://s2geometry.io/devguide/s2cell_hierarchy"]{S2 Cell
Hierarchy} to understand how geoids work.  What this library calls geoids, are
called cells by the S2 libary.

Also, while inspired by the S2 library, this Racket module is an independent
implementation that does not generate compatible IDs and does not aim to have
the same API and functionality as that library.

@itemlist[

@item{A geoid is a 64 bit integer between @racket[first-valid-geoid] and
      @racket[last-valid-geoid].  In particular, @racket[0] is not a valid
      geoid (they actualy start at @racket[1])}

@item{Geoids which are close together represent geographic locations which are
      close together, but the reverse is not true: there are geographic
      locations which are close together, but their geoids are very different}

@item{The projection method, while not uniform accross the globe, it only has
      a small distortion and there are no singularities anywhere, including at
      the poles}

@item{Geoids split the earth surface at different levels: the highest level is
      level @racket[30], where the earch surface is split into 6 faces.  At
      each lower level, the geoids are split into four, for example at level
      @racket[29], each of the 6 faces are split into four, producing 24 total
      geois.  The subdivision goes on until level @racket[0] is reached.  At
      level @racket[0], each geoid represents a patch of earth of
      approximately 0.7 cm@superscript{2}.}

]

@section{Converting to and from Latitude/Longitude}

@defproc*[([(lat-lng->geoid [lat real?] [lng real?] [#:level level exact-integer? 0]) exact-integer?]
          [(geoid->lat-lng [geoid exact-integer?]) (values real? real?)])]{

  Convert a geographic coordinate (latidude, longitude) to and from a geoid.
  @racket[lat-lng->geoid] will return a geoid at the specified @racket[level]
  of precision, by default at the highest level of precision, level
  @racket[0]. @racket[geoid->lat-lng] will accept geoids at any level of
  precision.

  @racket[geoid->lat-lng] will return the latitude/longitude coordinates
  corresponding to the center of the geoid and this will introduce an error in
  the conversion from latitude/longitude to geoid and back.  For geoids at
  level @racket[0], emprirical testing showed this to be less than 9.75
  millimeters, which is sufficient for the type of applications intended for
  this libary.  Note that this error is not constant accross the globe, and
  9.75 millimeters is the maximum error seen.

}

@defproc[(lat-lng-rect [geoid exact-integer?]) (values real? real? real? real?)]{

  Return the latitude / longitude rectangle which encloses the @racket[geoid].
  The four values returned are: minimum latitude, minimum longitude, maximum
  latitude, maximum longitude.  The bounds are slightly extended, to ensure
  all leaf geoids are inside the bounding box.

  The bounding box encloses the @racket[geoid] minimally, but geoids and
  bounding boxes don't overlap exactly, so the bounding box will always be
  bigger than then @racket[geoid] and there will be other geoids which are
  inside the bounding box, but not in the geoid.

}

@section{Storing and Retrieving from a SQLite database}

@defproc*[([(geoid->sqlite-integer [geoid exact-integer?]) exact-integer?]
          [(sqlite-integer->geoid [i exact-integer?]) exact-integer?])]{

  Convert a @racket[geoid] into an integer suitable for storage into a SQLite
  database, or convert an integer stored in the database back into a geoid.

  Geoids are unsigned 64 bit values, and more than half of them have the
  highest bit set.  SQLite will store numbers as signed 64 bits and convert
  unsigned 64 bit numbers greater than the maximum signed 64 bit value to
  floating point numbers, loosing precision.  This means that geoids cannot be
  stored directly into a SQLite database.

  These pair of functions subtract 2@superscript{63} from the geoid (or add
  that value back) to make sure the value is stored correctly.  The ordering
  of the geoids is preserved, so they can still be used to index geograpic
  information.

}

@section{Generating Test Data}

@defproc[(random-geoid [level exact-integer?] [#:parent parent (or/c exact-integer #f)]) exact-integer?]{

  Generate a random geoid at the specified @racket[level].  If @racket[parent]
  is specified, the level has to be lower (mode detailed) than the
  @racket[parent] and a geoid which is inside the @racket[parent] will be
  generated.

  This function is intended to generate geoids for unit tests.

}

@defproc[(leaf-outline [geoid exact-integer?]
                       [#:steps steps (and/c integer? positive?) 10]
                       [#:closed? closed? boolean? #t])
         (listof exact-integer?)]{

  Return a sequence of leaf geoids representing the outline of @racket[geoid].
  This can be used to obtain the outline of a geoid to display on a map.

  @racket[steps] specifies the number of geoids ot put on each side of the
  rectangle, while @racket[closed?] specifies if the first geoid should be
  duplicated as the last element in the list, to close the loop.

  As with @racket[leaf-corners], geoids are placed in counter-clockwise order,
  but there is no guarantee about the start corner of the loop.

}

@section{API Documentation}

@defproc*[([(first-valid-geoid) exact-integer?]
          [(last-valid-geoid) exact-integer?])]{

  Return the first and last valid integers that represent geoids.  All
  integers between these two numbers are valid geoids, but they can be at
  different levels.  Note that @racket[0] is not a valid geoid.

  See @racket[geoid-stride] to determine the amount to add to a geoid to
  obtain the next geoid at the same level.

}

@defproc[(sentinel-geoid) exact-integer?]{

  Returns an integer which is one bigger than the largest valid geoid (as
  returned by @racket[last-valid-geoid]).  The returned value is not a valid
  geoid.

  This can be used to create half-open geoid ranges, for example by
  @racket[leaf-span].

}

@defproc[(valid-geoid? (geoid exact-integer?)) boolean?]{

  Return true if @racket[geoid] is a valid geoid.  This is the same as testing
  that the @racket[geoid] is between @racket[first-valid-geoid] and
  @racket[last-valid-geoid], but shorter to type.

}

@defproc[(geoid-level [geoid exact-integer?]) (integer-in 0 30)]{

  Return the level of this geoid -- geoids can represent areas at increasing
  level of detail.  The highest resolution is at level @racket[0], and going
  up, geoids become 4 times bigger at each level.  Level @racket[30] is the
  top level, where the entire Earth surface is divided into @racket[6] faces.

  Geoids produced by @racket[lat-lng->geoid] are at level @racket[0] and you
  can use @racket[enclosing-geoid] to obtain a geoid at a higher level.

}

@defproc[(geoid-stride [geoid exact-integer?]) exact-integer?]{

  Return the integer step that must be added to a geoid to obtain another
  geoid at the same level.  This can be used to generate valid geoids in
  sequence.  Note that the geoids at the highest resolution (level @racket[2])
  have a stride of @racket[2], so you cannot simply incremenet geoids.

}

@defproc[(enclosing-geoid [geoid exact-integer?] [level (integer-in 0 30)])
         exact-integer?]{

  Return the geoid at @racket[level] that contains the specified
  @racket[geoid].  This function can be used to obtain the geoid at a lower
  resolution which contains a given point or geoid.

}

@defproc[(split-geoid [geoid exact-integer?])
         (list/c exact-integer? exact-integer? exact-integer? exact-integer?)]{

  Return the four geoids at the lower level into which the current
  @racket[geoid] can be split.  An error is raised, if the supplied geoid is a
  leaf geoid, at level @racket[0].

  The returned geoids are not in any particular order.

}

@defproc[(leaf-geoid? [geoid exact-integer?]) boolean?]{

  Return @racket[#t] if @racket[geoid] is a geoid at level 0 (highest level of
  detail).  This is equivalent to checking if @racket[geoid-level] returns
  @racket[0], but faster.

}

@defproc[(leaf-span [geoid exact-integer?]) (values exact-integer? exact-integer?)]{

  Returns the half-open geoid range that are valid geoids contained in
  @racket[geoid].  The first value is the smallest leaf geoid which is inside
  this @racket[geoid], the second value is either the smallest leaf geoid
  which is @bold{not} part of @racket[geoid], or the @racket[sentinel-geoid],
  whichever is smaller.

  All geoids which are inside this geoid, regardless of level, are contained
  in the returned number range, so this range can be used to check if any
  geoid is inside this one.

  The leaf span returned by this function can be used to search for geoids in
  an SQL query, however, if you do that, make sure you adjust them, as noted
  in the SQLite section above.

}

@defproc[(leaf-span* [geoids (list-of exact-integer?)]) (list-of (cons/c exact-integer? exact-integer?))]{

  Return a list of half open geoid ranges which define the valid geoids
  contained in all the geoids from @racket[geoids] list.  This is equivalent
  to calling @racket[leaf-span] for each individual geoid in @racket[geoids]
  and merging all the adjacent ranges together.

  This function can be used to create more efficient queries if a geoid is
  inside a list of geoids.  A common use case is to use
  @racket[adjacent-geoids] to obtain the neighbours of a geoid and using
  @racket[leaf-span*] to find the ranges for all geoids in this neighbourhood.

  }

@defproc[(contains-geoid? [this-geoid exact-integer?] [other-geoid exact-integer?]) boolean?]{

  Return true if the @racket[other-geoid] is geographically inside
  @racket[this-geoid].

  This a convenient function, but if you need to check lots of geoids, this
  will be slower than obtainging the @racket[leaf-span] of @racket[this-geoid]
  and checking of other geoids are inside the returned range.

}

@defproc[(leaf-corners [geoid exact-integer?])
         (list/c exact-integer? exact-integer? exact-integer? exact-integer?)]{

  Return the four leaf geoids which represent the corners of this
  @racket[geoid].  The corners are returned in couter-clockwise order, but
  there is no guarante of which one is first (i.e. there is no guarantee that
  the list will start with the top-left corner)

}

@defproc[(adjacent-geoids [geoid exact-integer?])
         (list-of exact-integer?)]{

  Return the adjacent geoids which border @racket[geoid].  The returned geoids
  will be at the same level as @racket[geoids].

  Normally, 8 geoids are returned, but only 7 are returned if @racket[geoid]
  is in a corner of a face and only 4 geoids if it is a face level geoid.

}

@defproc*[([(approximate-area-for-geoid [geoid exact-integer?]) real?]
           [(approximate-area-for-level [level (between/c 0 30)]) real?])]{

  Return the approximate area, in square meters, of @racket[geoid] or
  @racket[level].  The area is calculated by dividing the earth surface by the
  number of geoids at that level and it is approximate because there will be
  geoids with smaller area and larger area than this at each level.

  These functions can be used to get a general idea of the surface covered by
  a geoid.

}

@defproc[(distance-between-geoids [g1 exact-integer?] [g2 exact-integer?]) real?]{

  Return the distance, in meters, on the Earth surface between the center of
  geoids @racket[g1] and @racket[g2].  For leaf geoids this is a good
  appoximation for the distance between the locations represented by these
  geoids, since the size of a leaf geoid is approximately 8.5 millimeters.

}

@defproc[(distance-from-geoid [g exact-integer?]) (-> exact-integer? real?)]{

  Return a function which can be used to calculate the distance on the Earth
  surface between the geoid @racket[g] and another geoid.  The returned
  function will accept a single geoid as an argument and will return the
  distance between that geoid and @racket[g].

  If you need to calculate the distance between a single geoid and several
  others, it might be faster to use this function to construct a "distance
  function".

}

@section{Waypoint Alignment}
@defmodule[geoid/waypoint-alignment]

@defproc[(waypoint-alignment-cost [path1 (vectorof exact-integer?)] [path2 (vectorof exact-integer?)]) real?]{

  Return a number representing how similar @racket[path1] is to
  @racket[path2].  Both paths are vector of integers representing geoids.

  A smaller cost indicates that the two paths are more similar.  Ideally, the
  cost of a path against itself should be zero, but, due to floating point
  errors, this is a small positive number.

  This function returns the
  @hyperlink["https://en.wikipedia.org/wiki/Dynamic_time_warping"]{Dynamic
  Time Warping} cost of the two paths.

  @bold{NOTE}: the alignment cost will depend not only on how close the two
  paths are to each other, but also on the length of the paths, so it is up to
  the caller of the funtion to decide how to interpret the resulting cost and
  determine if the two paths are the same or not.

}

@section{Region Tiling}
@defmodule[geoid/tiling]

The @racket[geoid/tiling] module contains functions to calculate geoid
coverings for regions covering the earth surface.  A geoid covering is a list
of geoids which are inside a region.  This functionality can be used to
implement fast region containment tests for geographic data.

@defclass[region% object% ()]{

  Region objects define a region on the earth surface and are used for tiling,
  but otherwise they are opaque objects.

}

@defproc[(make-spherical-cap [lat real?] [lon real?] [radius real?])
         (is-a?/c region%)]{

  Create a @racket[region%] which is a circle around the
  @racket[lat]/@racket[lon] point, where @racket[radius] specifies the radius
  of the circle, in meters.  The circle extends around the sphere representing
  Earth, thus the name "spherical cap".

}

@defproc[(make-open-polyline [track (listof (vector real? real?))])
         (is-a?/c region%)]{

  Create a @racket[region%] representing the lines connecting lat/lon points
  in @racket[track].

  This is not strictly a "region", but sometimes it is useful to determine the
  geoid covering for a stretch of road, and this can be done by defining the
  road as a sequence of points and tiling it with geoids.

}

@defproc[(make-closed-polyline [track (listof (vector real? real?))]
                               [#:ccw? ccw? #t])
         (is-a?/c region%)]{

  Create a @racket[region%] representing the inside region defined by the
  lines connecting lat/lon points in @racket[track].

  On a sphere (which the geoid library uses as a model for Earth), each closed
  list of points defines two regions, since the entire sphere surface is
  finite.  The simplest way to visualize this is with a sequence of points
  along the Earths equator, which defines two regions: the northen hemisphere
  and the southern hemisphere.  The @racket[ccw?] parameter specifies which
  region is the "inside" region defined by the @racket[track]: when
  @racket[ccw?] is @racket[#t], the inside is the region which has the points
  in the counter-clockwise order, or to the left of the segments defined by
  the points as we go around the @racket[track].

  See also @racket[guess-winding-order] for determining the winding order of a
  sequence of points.

  The closed polyline defined by @racket[track] cannot have segments that
  intersect each other.

}

@defproc[(join-regions [r (is-a?/c region%)] ...+)
         (is-a?/c region%)]{

  Create a new region which is the union of the regions passed in as input.

}

@defproc[(intersect-regions [r (is-a?/c region%)] ...+)
         (is-a?/c region%)]{

  Create a new region which is the intersection of the regions passed in as
  input.

  This function can be used to define regions from GeoJSON "polygon with
  holes" inputs -- since the inner polygons have oposite winding order from
  the outer polygon, they can be created with the same winding order as the
  outer one, thus defining the opposite, or outer region, and intersecting
  them with the outer region.

}

@defproc[(subtract-regions [r1 (is-a?/c region%)] [r2 (is-a?/c region%)])
         (is-a?/c region%)]{

  Create a region containg the points which are inside @racket[r1] but not
  inside @racket[r2].

}

@defproc[(guess-winding-order [track (listof (vector real? real?))])
         (or/c 'cw 'ccw)]{

  Guess the winding order (clockwise or counter-clockwise) for the polygon
  defined by @racket[track], which is a list of latitude/longitude points. The
  guess works for "usual" polygons which are much smaller than a hemisphere
  and assumes that the "inside" is the smaller part.

  Regions deifned in GeoJSON objects don't specify a winding order, so this
  function can be used to determine the winding order of polygons inside
  GeoJSON objects, which can be passed on to @racket[make-closed-polyline].

}

@defproc[(geoid-tiling-for-region [r (is-a?/c region%)]
                                  [min-level integer?]
                                  [max-level integer?])
         (listof integer?)]{

  Return a list of geoids which cover the region @racket[r].  The region is
  covered with geoids of level @racket[max-level] first, than refined along
  the region borders with geoids down to @racket[min-level].

  See @racket[geoid-level] for a discution on the geoid levels.

  This function calls @racket[coarse-geoid-tiling-for-region] than
  @racket[refine-geoid-tiling-for-region] for each geoid which intersects the
  region, and returns all the geoids in a single list.

}

@defproc[(coarse-geoid-tiling-for-region [r (is-a?/c region%)]
                                         [level integer?])
         (values (listof integer?) (listof integer?))]{

  Return two lists of geoids at @racket[level] which cover the region
  @racket[r].  All geoids will be the same level.  The first list are all the
  geoids which are entirely inside the region @racket[r], while the second
  list contains all the geoids which intersect the region.

  See @racket[geoid-level] for a discution on the geoid levels.

}

@defproc[(refine-geoid-tiling-for-region [r (is-a?/c region%)]
                                         [geoid integer?]
                                         [min-level integer?])
         (listof integer?)]{

  Return a list of geoids which cover the region @racket[r] up to
  @racket[min-level] obtained by recursively splitting the input
  @racket[geoid].

  See @racket[geoid-level] for a discution on the geoid levels.

}

@section{Geodesy Calculations}
@defmodule[geoid/geodesy]

The @racket[geoid/geodesy] module contains functions to calculate distances
and bearings between points on the Earth surface.  Two calculation modes are
provided: one that approximates Earth as an ellipsoid, by default
@racket[wgs84], which is the model used by the GPS satelites, and another
which approximates the Earth as a sphere.

@defproc[(ellipsoid? [e any/c]) boolean?]{

  Return @racket[#t] if @racket[e] is an ellipsoid created by
  @racket[make-ellipsoid].

}

@defproc[(make-ellipsoid [major real?] [minor real?]) ellipsoid?]{

  Create an ellipsoid with the @racket[major] and @racket[minor] semi-axes.
  Ellipsoids can be used as values for @racket[geodesy-ellipsoid] or to
  individual functions in this module.  The @racket[wgs84] ellipsoid, used by
  the GPS satelites, is already defined in this module and it is the default.

}

@defthing[wgs84 ellipsoid?]{

  The @hyperlink["https://en.wikipedia.org/wiki/World_Geodetic_System"]{WGS84}
  ellipsoid.

}

@defparam[geodesy-angle-mode mode (or/c 'degrees 'radians) #:value 'degrees]{

  A parameter that specifies the type of angles which are passed to the
  functions in this module as well as the type of angle which are returned for
  bearings.  Valid values are @racket['degees], which means that latitude,
  longitudes are in degrees or @racket['radians].

}

@defparam[geodesy-ellipsoid ellipsoid (or/c #f ellipsoid?) #:value wgs84]{

  A parameter that specifies the ellipsoid to use for the calculations by the
  functions in this module.  By default, the @racket[wgs84] ellipsoid is used
  and a value of @racket[#f] means that the calculations are done assuming a
  sperical Earth model.

}

@defproc[(distance-between [lat1 real?] [lon1 real?]
                           [lat2 real?] [lon2 real?]
                           [#:angle-mode m (or/c 'degrees 'radians) (geodesy-angle-mode)]
                           [#:ellipsoid e (or/c #f ellipsoid?) (geodesy-ellipsoid)])
         real?]{

  Return the distance between two points on Earth along the great circle,
  which is the shortest distance.  The points are identified by latitude and
  longitude, which are angles.

  @racket[m] specifies if the latitude and longitude angles are specified in
  degrees or radians, while @racket[e] specifies the ellipsoid to use for the
  calculation. If @racket[e] is @racket[#f], the calculation is done assuming
  Earth is a sphere.

}

@defproc[(initial-bearing [lat1 real?] [lon1 real?]
                          [lat2 real?] [lon2 real?]
                          [#:angle-mode m (or/c 'degrees 'radians) (geodesy-angle-mode)]
                          [#:ellipsoid e (or/c #f ellipsoid?) (geodesy-ellipsoid)])
         real?]{

  Return the initial bearing for traveling between two points on Earth along
  the great circle.  The points are identified by latitude and longitude,
  which are angles and the bearning is also an angle, where @racket[0] is the
  North direction and the angle moves clockwise.

  Note that, when traveling along a great circle the bearing will not remain
  constant.

  The @racket[m] and @racket[e] parameters are the same as for
  @racket[distance-between].

}

@defproc[(final-bearing [lat1 real?] [lon1 real?]
                        [lat2 real?] [lon2 real?]
                        [#:angle-mode m (or/c 'degrees 'radians) (geodesy-angle-mode)]
                        [#:ellipsoid e (or/c #f ellipsoid?) (geodesy-ellipsoid)])
         real?]{

  Similar to @racket[initial-bearing], but return the final bearing for
  traveling between two points on Earth along the great circle.

}

@defproc[(destination-point [lat real?] [lon real?]
                            [bearing real?] [distance real?]
                            [#:angle-mode m (or/c 'degrees 'radians) (geodesy-angle-mode)]
                            [#:ellipsoid e (or/c #f ellipsoid?) (geodesy-ellipsoid)])
         (values real? real?)]{

  Return the destination point for traveling from the point @racket[lat],
  @racket[lon] along a great circle with an initial @racket[bearing] for a
  given @racket[distance].

  The @racket[m] and @racket[e] parameters are the same as for
  @racket[distance-between].

}

@defproc[(midway-point [lat1 real?] [lon1 real?]
                       [lat2 real?] [lon2 real?]
                       [#:angle-mode m (or/c 'degrees 'radians) (geodesy-angle-mode)]
                       [#:ellipsoid e (or/c #f ellipsoid?) (geodesy-ellipsoid)])
         (values real? real?)]{

  Return the latitude/longitude point that it half-way between the points
  @racket[lat1], @racket[lon1] and @racket[lat2], @racket[lon2], along the
  great circle arc, that is along the shorted distance path between the two
  points.

  The @racket[m] and @racket[e] parameters are the same as for
  @racket[distance-between].

}
