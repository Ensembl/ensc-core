/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __ERROR_H__
#define __ERROR_H__

void     Error_trace(char *Routine, char *Info);
#ifdef __hpux
  void     Error_write();
#else
  void     Error_write(int errorno,char *routine,int level,char *extra,...);
#endif

/* Error levels */
#define ERR_MILD        1
#define ERR_SEVERE      2
#define ERR_FATAL       3

/* IO errors */
#define EFOPEN 		1
#define EFGETSNULL	2
#define EFGETSZERO	3
#define EGETTOKZERO     4
#define EMAXLENTOK      5
#define EFTELL          6
#define EEOF            7
#define ETEMPNAM        8
#define WEMPTYFILE      9 
#define EOPENDIR        10

/* Pointer Errors */
#define ENULLPOINTER	11
#define EMALLERR	12
#define ECALLERR	13
#define EVECERR         14
#define EREALLERR       15

/* AA errors */
#define EUNKAA          16

/* String Errors */
#define EDIFFSTRING     18
#define ESTRNCATNULL    19
#define EZEROSTRING	20
#define ESSCANFNUM	21
#define ESTRTOLONG	22
#define EZEROLEN        23
#define ESTRNCPYNULL    24
#define ESTRTOSHORT     25
#define ESSCANF         26
#define EVSPRINTF       27
#define ESTRTOL         28

#define ESINGMAT        29
/* Array errors */
#define EZEROARRAY      30

/* Fitting errors */
#define EZEROMAT        31
#define EROTSPECIAL     32
#define EINVALNAT       33
#define EMAXITER        34

/* Application errors */
#define EUSAGE          35
#define ENOPDB          36
#define ENUMSEQ         37
#define EORDERERR       38
#define EINVALCUT       39
#define EPDBINCONSIST   40
#define ENOORDERRES     41
#define EXTRAAT         42
#define EMISSAT         43
/* SMJS #define ENOATOM         44 (redefined later) */
#define EZEROSEQ        45
#define EMISMATCH       46
#define ENODEL          47
#define EXTRARES        48
/* SMJS #define ENOSEQ          49 (redefined later) */
#define ENUMTOOLONG     50
#define EINCLEN         51
#define EINCZERO        52
#define ENORES          53
#define ENOCOFGAT       54
#define EDIFFNUMDIFF    55
#define EZEROLENSEQ     56
#define ENORANGE        57
#define EPDBRESMISSING  58
#define ERANGELIMITS    59
#define EINCERROR       60
#define EMATFILEERR     61
#define ECHNNOTFOUND    62
#define ESEQPOSRANGE    63
#define ESEQWINGAP      64
#define EDIFFLENSEQ     65
#define EZERORES        66
#define EINSUFRES       67

/* Parse Defines */
#define ERRNOUNIQ       68
#define ERRNOSTR        69
#define ERRZEROLENSTR   70
#define ERECLEV         71
#define EBLKSTRT        72
#define ENOARG          73
#define EPARSEFF        74
#define ENOBRACE        75
#define ENOKEYVAL       76 // Note changed from ENOKEY to avoid conflict with errno.h
#define EKEYID          77
#define EKEYDEF         78
#define EARG            79
#define ELISTLIM        80
#define EUNKTYPE        81
#define EARGDEF         82
#define EARGTYPE        83
#define ETYPEFUNC       84
#define EARGLISTNUM     85
#define EINT            86
#define EFLOAT          87
#define ELISTFORMAT     88
#define EDEFFLAG        89
#define EUNKSTYPE       90
#define EFALLTHRU       91
#define EMODE           92
#define ESOURCE         93
#define ESYNTAX         94
#define ESTRNOTFND      95
#define EINVSTYPE       96
#define ENOCOMMAND      97
#define ENOTIMPL        98
#define EARGRANGE       99
#define EMAXARRAY      100
#define EAMBIGKEY      101
#define EKEYNOTFND     102
#define EKEYNAME       103
#define EHARDARG       104
#define ENOREQARG      105
#define ENOARGVAL      106
#define ENODEFAULT     107
#define ENUMRANGE      108
#define WNOHIST        109
#define ENUMARG        110
#define ERANGEFORMAT   111
#define ECONTINUE      112
#define EQUOTE         113
#define WEXTRANEOUS    114
#define EMISSITEM      115
#define WMISSCHAR      116
#define ENUMTOK        117
#define EDUPLICATE     118
#define EMISSCHAR      119
#define EARGVALUE      120

/* wild card errors */
#define ESTRTBRACE     121
#define EENDBRACE      122
#define ELINKCHAR      123
#define EMISSCOMMA     124
#define EBRACE         125

/* Inconsistancies */
#define EINCONSIST     129
#define ENOMATCH       130
#define EACCVIOL       131

/* List errors */
#define ENOENTRY       132
#define ENOOBJ         133
#define ETOOMANYOBJ    134
#define EOBJTYPE       135
#define E2WAYLIST      136
#define ENONUNIQUENAME 137

/* Sequence errors */
#define ENOLINE        140
#define EINSURFSEQ     141
#define EPOSVAL        142
#define ESEQFORMAT     143
#define ENOCHAIN       144
#define ESEQNUM        145
#define ENOSEQ         146
#define EDIFFNUMCH     147
#define EDIFFNUMSEQ    148
#define EFREQTOT       149


/* Selection errors */
#define ESELFORMAT     150
#define ENUMFORMAT     151
#define ENUMOBJ        152
#define ENUMATOM       153
#define EDIFFNUMAT     154
#define ENOSELECT      155
#define EDIFFNUMRES    156
#define EDIFFCHAIN     157
#define ESELNUM        158
#define ENOSELAT       159
#define EOPERATOR      160
#define ECMPSELFMT     161
#define ENOFREEGRP     162
#define ENSELRES       163
#define EDIFFNUMOBJ    164

/* Math errors */
#define ERANDNUM       170

/* Config errors */
#define ECONFIGFORMAT  180
#define ECONFSEQ       181
#define ECONFPDB       182

/* Stream errors */
#define EUNKNOWNSTREAM 190
#define EINVALIDSTREAM 191

/* Atom and Residue Errors */
#define EUNKATOM       200
#define EMISSINGAT     201
#define ENOATOM        202
#define EMISSINGRES    203
#define ECOUNTRES      204

/* Torsion Errors */
#define ENOSUCHTOR     210
#define EEQUIVAT       211
#define ENOMORERES     212

/* File Errors */
#define ESAMEFILE      220
#define EINVMODE       221
#define EFFILEFMT      222

/* Discover I/0 Errors */
#define ENOTSPACE      230
#define EBLOCKLEN      231
#define EFFFORMAT      232
#define EZEROSECTION   233
#define EFFTYPE        234
#define ECARFORMAT     235
#define EMDFFORMAT     236
#define ERESLIBFMT     237

/* Discover molecule errors */
#define EUNKMOLETYP    240
#define ENOMOL         241
#define EBONDMISSING   242
#define EBONDORDER     243
#define ENOPOT         244
#define ENOBOND        245
#define EINVALIND      246
#define ENOANGLE       247
#define ENOOOP         248
#define ENOTORSION     249

#define EINVALOOPAT    250
#define EMISSINGANGLE  251
#define ENOOOPOOP      252
#define ENOANGLEANGLE  253
#define ENUMPOT        254
#define ENOSWITCHATOM  255
#define EZERODIST      256
#define WNOATOMTYPE    257
#define EFILETYPE      258
#define EBONDNOBJ      259

/* residue library errors */
#define EBONDNATOM     260
#define EGENRES        261
#define EEMBEDTER      262
#define EDUMMYATOM     263

/* accessibility errors */
#define ENUMVERTEX     270
#define ENOICOSA       271

/* Kabat format errors */
#define ERESTYPE       280
#define EKABFORMAT     281
#define EKABUNSUP      282
#define EKABFFFORMAT   283
#define EKABINSERT     284
#define WNOTKABAT      285
#define ENOTKABAT      286

/* command argument handler errors */
#define WOPTCHAR       290
#define EOPTFORMAT     291
#define EOPTUNK        292
#define EOPTNUM        293
#define EOPTTYPE       294
#define WOPTMISS       295
#define EOPTNOTFOUND   296
#define EOPTERR        297

/* More system errors */
#define EGETENV        300
#define ESTRTOD        301
#define ESYSTEM        302
#define EEXISTS        303
#define WEXISTS        304
#define EFPRINTF       305
#define EFPUTS         306
#define EGETCWD        307
#define EREALPATH      308
#define ECHDIR         309

/* CGA conformation errors */
#define ECGAPDB        310
#define ECGACAR        311
#define ECGACGA        312
#define EINSURFCONF    313
#define ECGANCONS      314
#define EDIFFATTYPE    315
#define WINSURFCONF    316
#define ECGANOCONF     317
#define ENOCGAPDB      318
#define ENONSEQCGA     319

/* torsion errors */
#define ETORTYPE       320
#define ERINGTOR       321
#define EMISSANCEST    322

/* Conformation generation errors */
#define EANGFORMAT     330
#define ETORRESIDUE    331
#define EWRONGOBJ      332
#define ENSELATOM      333
#define EATOMTYPE      334
#define EMISSINGDOF    335
#define ENOPREVRES     336
#define ENODERIV       337
#define EINVALITER     338
#define ENOENERGY      339

/* Graphics errors */
#define EGRNOTINIT     350

/* Number file errors */
#define ENSCHEMEFORMAT 360
#define ESETSCHEME     361
#define EWRONGNUMTYPE  362

/* More Sequence errors */
#define ETENBFORMAT    370
#define ELENERR        371
#define EUNKSPEC       372
#define EPROMPTCHN     373
#define EEMBLFORMAT    374
#define EWRONGNUMCHN   375
#define EVARINUSE      376
#define EDATATYPE      377
#define EMISSCOMMENT   378

/* EcoString and CHash errors */
#define EECOSTR        380
#define ECHASH         381

/* More Discover errors */
#define EINCFFIELD     390

/* OxBench errors */
#define ENOXBFAM       400
#define EOXBFLAGS      401
#define ESTAMPFMT      402
#define ESPLITOXBNAME  403

/* PBuffer errors */
#define EPBEMPTY       410
#define EPBFULL        411
#define EPBID          412

/* More system errors */
#define WGETENV        420

/* MySQL errors */
#define EMYSQLCONN     430

/* OBDA errors */
#define EOBDA          440

/* !!You must remember to add the error message to Message.h and increase the */
/*   total number of errors set in Message.h */

#ifndef __MESSAGE_H__
extern char *message[];
#endif /* __MESSAGE_H__ */

#ifndef NOEXTERN
extern int ExitLevel;
extern int LastErrNo;
extern int LastELev;
#endif

#endif /* __ERROR_H__ */
