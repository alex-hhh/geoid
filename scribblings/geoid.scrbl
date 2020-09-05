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

@section{Storing geoids in a SQLite database}

Geoids use all 64 bits, and more than half of them have the highest bit set.
SQLite will store numbers as signed 64 bits and convert unsigned 64 bit
numbers to 64 bit floating points.

To store geoids in a SQLLite database, you need to subtract 2@superscript{63}
from it, to convert it to a signed number.  This will preserve the ordering
properties of the geoid.  Just remember to add back 2@superscript{63} back to
it before using a geoid retrived from the database this way.

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

@defproc*[([(lat-lng->geoid [lat real?] [lng real?]) exact-integer?]
          [(geoid->lat-lng [geoid exact-integer?]) (values real? real?)])]{

  Convert a geographic coordinate (latidude, longitude) to and from a geoid.
  The returned geoid will be at the highest level of precision (level
  @racket[0]).

  The conversion from latitude/longitude to geoid and back has a small amount
  of error, emprirical testing showed this to be less than 7 millimeters,
  which is sufficient for the type of applications intended for this libary.
  Note that this error is not constant accross the globe, and 7 millimeters is
  the maximum error seen.

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
