#include "BaseComparaAdaptor.h"

void BaseComparaAdaptor_init(BaseComparaAdaptor *bca, ComparaDBAdaptor *dba, int adaptorType) {
  bca->dba = dba;
  bca->adaptorType = adaptorType;
  bca->prepare = BaseAdaptor_prepare;
}
