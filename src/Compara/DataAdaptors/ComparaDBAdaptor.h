#ifndef __COMPARADBADAPTOR_H__
#define __COMPARADBADAPTOR_H__

#include "BaseDBAdaptor.h"
#include "DBConnection.h"
#include "ComparaAdaptorTypes.h"
#include "EnsC.h"
#include "StringHash.h"

struct ComparaDBAdaptorStruct {
  BASEDBADAPTOR_DATA
  StringHash *genomes;
};

ComparaDBAdaptor *ComparaDBAdaptor_new(char *host, char *user, char *pass, char *dbname,
                                       unsigned int port, char *confFile);

#endif
