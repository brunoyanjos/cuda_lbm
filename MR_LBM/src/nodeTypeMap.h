#ifndef __NODE_TYPE_MAP_H
#define __NODE_TYPE_MAP_H

#include <builtin_types.h>
#include <stdint.h>

// DIRECTION DEFINES 00000000

#define SOLID_NODE (0b00000000)

#define NORTH (0b00000011)
#define SOUTH (0b00001100)
#define WEST (0b00001010)
#define EAST (0b00000101)

#define NORTH_WEST (0b00000010)
#define NORTH_EAST (0b00000001)
#define SOUTH_WEST (0b00001000)
#define SOUTH_EAST (0b00000100)

#define INNER_BOUNDARY (0b00010000)

#define BULK (0b00001111)

#endif // !__NODE_TYPE_MAP_H
