#include "RepeatConsensusAdaptor.h"
#include "BaseAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"


RepeatConsensusAdaptor *RepeatConsensusAdaptor_new(DBAdaptor *dba) {
  RepeatConsensusAdaptor *rca;

  if ((rca = (RepeatConsensusAdaptor *)calloc(1,sizeof(RepeatConsensusAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for RepeatConsensusAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)rca, dba, REPEATCONSENSUS_ADAPTOR);

  return rca;
}

RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *rca, int64 dbID) {
  RepeatConsensus *repeatConsensus;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if (IDHash_contains(ca->chrCache,dbID)) {

    repeatConsensus = IDHash_getValue(ca->chrCache, dbID);

  } else {
    sprintf(qStr,"SELECT repeat_consensus_id, name, length"
      " FROM repeat_consensus"
      " WHERE  repeat_consensus_id = "
      INT64FMTSTR, dbID);
  
    sth = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));
    sth->execute(sth);
  
    row = sth->fetchRow(sth);
    if( row == NULL ) {
      sth->finish(sth);
      return NULL;
    }
  
    repeatConsensus = RepeatConsensusAdaptor_repeatConsensusFromRow(ca, row);
    sth->finish(sth);
  }

  return repeatConsensus;
}
