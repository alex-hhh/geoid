#lang scribble/manual
@require[@for-label[geoid
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

@defproc*[([(lat-lng->geoid [lat real?] [lng real?]) exact-integer?]
          [(geoid->lat-lng [geoid exact-integer?]) (values real? real?)])]{

  Convert a geographic coordinate (latidude, longitude) to and from a geoid.
  The returned geoid will be at the highest level of precision (level
  @racket[0]), but @racket[geoid->lat-lng] will also accept geoids at a lower
  precision (a level higher than @racket[0]).

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
         (list-of integer?)]{

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


