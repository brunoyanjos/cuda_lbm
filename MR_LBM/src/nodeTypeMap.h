#ifndef __NODE_TYPE_MAP_H
#define __NODE_TYPE_MAP_H

#include <builtin_types.h>
#include <stdint.h>

// DIRECTION DEFINES 00000000

#define BULK (15)
// FACE
#define NORTH (3)
#define SOUTH (12)
#define WEST (10)
#define EAST (5) 
// CORNER
#define NORTH_WEST (2)
#define NORTH_EAST (1)
#define SOUTH_WEST (8)
#define SOUTH_EAST (4)
// IMMERSED

#define SOLID_NODE (0)

#define MISSING_DEFINITION (0b11111111111111111111111111111111)

#define DIRECTION_BITS (0b11111 << DIRECTION_OFFSET)

#endif // !__NODE_TYPE_MAP_H
