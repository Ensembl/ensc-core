/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __MESSAGE_H__
#define __MESSAGE_H__

/* General defines */
#define NUMERROR       441


char *message[NUMERROR] = {
	"Error: No error (check the message number)!! ",
	"Error: Unable to open file ",
	"Error: fgets() returned a NULL pointer ",
	"Error: Attempt to get a zero length string using fgets() ",
	"Error: Gettok failed to find a token",
	"Error: Maximum length of token exceeded",
	"Error: ftell() returned an error",
	"Error: End of file",
	"Error: Error creating temporary name",
	"Warning: Empty file",
/* 10 */
	"Error: Opening directory",
	"Error: Null pointer ",
	"Error: Failed in malloc. Variable name ",
	"Error: Failed in calloc. Variable name ",
	"Error: Failed allocating vector called from",
	"Error: Failed in realloc. Variable name ",
	"Error: Unknown amino acid ",
	"",
	"Error: Strings are different",
	"Error: strncat returned a NULL",
/* 20 */
	"Error: Zero length string ",
	"Error: sscanf() returned an incorrect number of arguments ",
	"Error: String too long. ",
	"Error: zero length ",
	"Error: strncpy returned a NULL",
	"Error: string too short",
	"Error: sscanf error",
	"Error: vsprintf error",
	"Error: Strtol didn't find an int",
	"Error: Singular matrix",
/* 30 */
	"Error: No elements in array",
	"Error: Singular matrix",
	"Error: Special case in fitting routine. No fitting done",
	"Error: Invalid number of atoms",
	"Error: Maximum number of iterations exceeded",
/* 35 */
	"Usage ",
	"Error: PDB file not read correctly",
	"Warning : Number of sequences read greater than 2",
	"Error: Error reading order file",
	"Error: Invalid cutoff ",
/* 40 */
	"Error: Inconsistancy ",
	"Error: Residue not found in order structure",
	"Warning: extra atoms found in residue",
	"Error: Missing atoms in residue while ordering",
	"Error: Atom not found for residue",
	"Error: No sequences",
	"Error: Sequence mismatch",
	"Error: Failed to create DEL residue ",
	"Error: Extra residues after aligned region",
	"Error: No sequence in file ",
/* 50 */
	"Error: String representing residue number is too long",
	"Error: String representing residue number is too long",
	"Error: String is zero long",
	"Error: No residues ",
	"Error: No atoms to calculate CofG for ",
	"Error: Different number of differences",
	"Error: Zero length sequence",
	"Error: No ranges",
	"Error: Failed to find PDB residue number",
	"Error: Range limits are wrong",
/* 60 */
	"Error: Error with PDB numbers",
	"Error: Reading MDM matrix file",
	"Error: Chain not found for a residue number",
	"Error: Sequence position out of range",
	"Warning: Sequence window contained a gap",
	"Error: Sequences are different lengths",
	"Error: Zero residues",
	"Error: Insufficient residues in sequence",
	"Error: Nonunique key",
	"Error: No strings",
/* 70 */
	"Error: String length zero",
	"Error: Recursion level error",
	"Error: Block start not found",
	"Error: No arguments for command",
	"Error: Parse file format error",
	"Error: No brace at start of block",
	"Error: No key at start of block",
	"Error: key Id error",
	"Error: Key default error",
	"Error: Argument should begin with OPTARG or REQARG or HARDARG\n",
/* 80 */
	"Error: Min and Max for type not correct",
	"Error: Unknown type",
	"Error: Argument default error",
	"Error: Argument not correct type",
	"Error: Type function is not correct for type",
	"Error: Number of arguments in list invalid",
	"Error: String is not an int",
	"Error: String is not an float",
	"Error: List format wrong",
	"Error: Default flag not NEVER or FIXED or LAST",
/* 90 */
	"Error: Unknown stream type for file",
	"Error: Fell through! ",
	"Error: Unknown stream mode",
	"Error: In source statement",
	"Error: Syntax error",
	"Error: String not found",
	"Error: Invald STREAM type for this action",
	"Error: No commands present",
	"Error: Command not implemented",
	"Error: Argument number not valid",
/* 100 */
	"Error: Maximum size of array exceeded",
	"Error: Ambiguous key",
	"Error: Key not found",
	"Error: Key name error",
	"Error: HARDARG requires an argument name",
	"Error: Extraneous tokens in command",
	"Error: Argument value not specified",
	"Error: No default argument",
	"Error: Number out of range",
	"Warning: There were no history commands to list",
/* 110 */
	"Error: Wrong number of arguments for command",
	"Error: Format of range string wrong",
	"Error: Continuation error",
	"Error: Quoted string error",
	"Warning: Extraneous characters",
	"Error: Missing list item",
	"Warning: Missing character",
	"Error: Wrong number of tokens",
	"Error: Duplicate entry",
	"Error: Missing character",
/* 120 */
	"Error: Argument value is not one of the allowed values",
	"Error: No string doesn't start with a brace",
	"Error: String didn't contain an end brace",
	"Error: Unknown link character",
	"Error: Missing comma",
	"Error: Missing brace",
	"",
	"",
	"",
	"Error: Inconsistancy",
/* 130 */
	"Error: No match",
	"Error: Outside array bounds",
	"Error: No entries",
	"Error: No objects match name",
	"Error: Too many objects",
	"Error: Unknown object type",
	"Error: Error in list",
	"Error: Non unique object name",
	"",
	"",
/* 140 */
	"Error: Couldn't find appropriate line in file",
	"Error: Insufficient sequences in file",
	"Error: Should be positive values",
	"Error: Sequence file format error",
	"Error: No such chain",
	"Error: No such sequence number in list",
	"Error: No sequences",
	"Error: Different numbers of chains",
	"Error: Different numbers of sequences",
	"Error: Frequencies do not add up to 1",
/* 150 */
	"Error: Selection string format (obj:chainname:resname:atom) wrong",
	"Error: Selector number format wrong",
	"Error: Number of objects wrong",
	"Error: Number of selected atoms wrong",
	"Error: Different numbers of atoms in selected objects",
	"Error: No Selection strings",
	"Error: Different numbers of residues",
	"Error: Different chains",
	"Error: Selection number out of range",
	"Error: No atoms selected",
/* 160 */
	"Error: Badly positioned operator",
	"Error: Complex selection format",
	"Error: No free selection groups",
	"Error: Number of selected residues",
	"Error: Different numbers of objects selected",
	"",
	"",
	"",
	"",
	"",
/* 170 */
	"Error: Random number generation problem",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 180 */
	"Error: Configuration file format error",
	"Error: Sequence and Configuration do not match",
	"Error: PDB and Configuration do not match",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 190 */
	"Error: Unknown stream type",
	"Error: Invalid file for stream",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 200 */
	"Error: Missing atom",
	"Error: Unknown atom type",
	"Error: No atom data",
	"Error: Missing residue",
	"Error: Counting residues between atoms",
	"",
	"",
	"",
	"",
	"",
/* 210 */
	"Error: No such torsion",
	"Error: No atom in equivalence table",
	"Error: No more residues",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 220 */
	"Error: Same file name",
	"Error: Invalid file mode",
	"Error: FFile format error",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 230 */
	"Warning: Non space character",
	"Error: Direct access block size error",
	"Error: Forcefield file format error",
	"Error: No non comment lines in section",
	"Error: File identification string incorrect",
	"Error: Car file format error",
	"Error: Mdf file format error",
	"Error: ResLib file format error",
	"",
	"",
/* 240 */
	"Error: Unknown molecule type",
	"Error: No molecules",
	"Error: Bond missing",
	"Error: Bond order error",
	"Error: No potential",
	"Warning: No bonds",
	"Error: Invalid index",
	"Warning: No angles",
	"Warning: No out_of_planes",
	"Warning: No torsions",

/* 250 */
	"Error: Invalid number of atoms for out_of_plane",
	"Error: Missing angle",
	"Warning: No out_of_plane cross terms",
	"Warning: No angle-angle cross terms",
	"Error: Wrong number of potentials",
	"Error: No switching atoms",
	"Error: Distance zero",
	"Warning: Atomtype not found",
	"Error: File type error",
	"Error: Inter molecule bond not allowed",

/* 260 */
	"Error: Only two atoms can make a bond",
	"Error: Generating coordinates for residue",
	"Error: Embedded terminal residue",
	"Error: Problem with dummy atoms",
	"",
	"",
	"",
	"",
	"",
	"",

/* 270 */
	"Error: Number of vertices wrong",
	"Error: No icosahedron",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 280 */
	"Error: Invalid Residue type",
	"Error: Kabat sequence format error",
	"Error: Unsupported format",
	"Error: Kabat file of files format error",
	"Error: Kabat inserted residue (AAIN) error",
	"Warning: Entry did not contain Kabat data",
	"Error: Entry is not Kabat data",
	"",
	"",
	"",
/* 290 */
	"Warning: option expects a char argument",
	"Error: Format of option string wrong",
	"Error: Unknown command line argument",
	"Error: Tried to get more options than command line arguments",
	"Error: Unknown format option",
	"Warning: Suspected missing argument before",
	"Error: option not found in options array",
	"Error: getarg messed up",
	"",
	"",
/* 300 */
	"Error: Unable to find enviroment variable",
	"Error: strtod() failed to decode a double",
	"Error: system() returned failure",
	"Error: File does not exist",
	"Warning: File does not appear to exist",
	"Error: fprintf() returned failure",
	"Error: fputs() returned failure",
	"Error: getcwd() failed",
	"Error: realpath() failed",
	"Error: chdir() failed",
/* 310 */
	"Error: CGA and PDB are inconsistent",
	"Error: CGA and CAR are inconsistent",
	"Error: CGA file and CGA structure are inconsistent",
	"Error: Insufficient CGA conformations",
	"Error: Number of constructed atoms different",
	"Error: Different atom types",
	"Warning: Insufficient CGA conformations",
	"Error: No CGA conformations",
	"Error: PDB for CGA is required but not set",
	"Error: Non sequential atoms in CGA",
/* 320 */
	"Error: Invalid torsion type string",
	"Error: Missing ancestors",
	"Error: Ring atoms selected for torsion",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 330 */
	"Error: Angle file format error",
	"Error: Didn't find residue",
	"Error: Wrong object selected",
	"Error: Wrong number of atoms selected",
	"Error: Wrong atom type",
	"Error: Missing DOF",
	"Error: No previous residues",
	"Error: No derivative data",
	"Error: Invalid iteration value",
	"Error: No energy data",
/* 340 */
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 350 */
	"Error: Graphics not initialised",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 360 */
	"Error: Number Scheme file format error",
	"Error: Set Scheme error",
	"Error: Wrong/No numbering set",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 370 */
	"Error: Ten block type sequence format error",
	"Error: Lengths error",
	"Error: Unknown type specifier",
	"Error: Choosing chain number",
	"Error: EMBL file format error",
	"Error: Wrong number of chains",
	"Error: Variable already in use",
	"Error: Wrong type of data",
	"Error: Comment missing for entry",
	"",
/* 380 */
	"Error: ECOSTRING handling error",
	"Error: CHASH handling error",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 390 */
	"Error: Incompatible forcefield object",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 400 */
	"Error: Number of OxBench sequences wrong",
	"Error: OXBench flags are incorrect",
	"Error: STAMP file format error",
	"Error: Incorrect format for OxBench split family name",
	"",
	"",
	"",
	"",
	"",
	"",
/* 410 */
	"Error: PBuffer empty",
	"Error: PBuffer full",
	"Error: PBuffer ID error",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 420 */
	"Warning: Environment variable not set",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 430 */
        "Error: Failed connecting to mysql database",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	"",
/* 440 */
        "Error: OBDA error"
};
#endif /* __MESSAGE_H__ */
