#ifndef __MAPPERDEFS_H__
#define __MAPPERDEFS_H__

#define MAPPER_FROM_IND 0
#define MAPPER_TO_IND   1

typedef enum CoordSystemEnum {
  ASSEMBLY_COORDS,
  RAWCONTIG_COORDS,
  CDNA_COORDS,
  GENOMIC_COORDS
} CoordSystemType;

#endif
