/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         scotchyyparse
#define yylex           scotchyylex
#define yyerror         scotchyyerror
#define yydebug         scotchyydebug
#define yynerrs         scotchyynerrs
#define yylval          scotchyylval
#define yychar          scotchyychar

/* First part of user prologue.  */
#line 1 "parser_yy.y"

/* Copyright 2004,2007,2008,2011,2014,2018 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : parser_yy.y                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the syntactic parser     **/
/**                which processes strategy strings.       **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 nov 1995     **/
/**                                 to     13 jun 1996     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to     27 feb 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 20 dec 2001     **/
/**                                 to     11 jun 2004     **/
/**                # Version 5.1  : from : 30 oct 2007     **/
/**                                 to     24 jul 2011     **/
/**                # Version 6.0  : from : 30 sep 2014     **/
/**                                 to     27 apr 2018     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define PARSER_YY

#include "module.h"
#include "common.h"

#undef INTEGER                                    /* In case someone defined them */
#undef DOUBLE

#include "parser.h"
#include "parser_ll.h"
#include "parser_yy.h"

/* #define SCOTCH_DEBUG_PARSER3 */
#ifdef SCOTCH_DEBUG_PARSER3
extern int                  yydebug;
#define YYDEBUG                     1
#endif /* SCOTCH_DEBUG_PARSER3 */

/*
**  The static and global definitions.
**  See also at the end of this file.
*/

static const StratTab *     parserstrattab;       /* Pointer to parsing tables          */
static Strat *              parserstratcurr = NULL; /* Pointer to current strategy node */
static StratParamTab *      parserparamcurr = NULL; /* Pointer to current parameter     */

extern unsigned int         parsermethtokentab[]; /* Pre-definition for stupid compilers */


#line 169 "y.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif


/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int scotchyydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    METHODNAME = 258,              /* METHODNAME  */
    PARAMNAME = 259,               /* PARAMNAME  */
    VALCASE = 260,                 /* VALCASE  */
    VALDOUBLE = 261,               /* VALDOUBLE  */
    VALINT = 262,                  /* VALINT  */
    VALSTRING = 263,               /* VALSTRING  */
    VALSTRAT = 264,                /* VALSTRAT  */
    VALPARAM = 265,                /* VALPARAM  */
    VALTEST = 266                  /* VALTEST  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif
/* Token kinds.  */
#define YYEMPTY -2
#define YYEOF 0
#define YYerror 256
#define YYUNDEF 257
#define METHODNAME 258
#define PARAMNAME 259
#define VALCASE 260
#define VALDOUBLE 261
#define VALINT 262
#define VALSTRING 263
#define VALSTRAT 264
#define VALPARAM 265
#define VALTEST 266

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 92 "parser_yy.y"

  char                      CASEVAL;              /* Case value          */
  StratTest *               TEST;                 /* Test type           */
  StratTestType             TESTOP;               /* Relational type     */
  double                    DOUBLE;               /* Double-precision    */
  INT                       INTEGER;              /* Integer             */
  char                      STRING[PARSERSTRINGLEN]; /* Character string */
  struct {
    const StratTab *        tabl;                 /* Current tables    */
    Strat *                 strat;                /* Current method    */
    StratParamTab *         param;                /* Current parameter */
  } SAVE;                                         /* Parameter type    */
  Strat *                   STRAT;                /* Strategy tree     */

#line 256 "y.tab.c"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE scotchyylval;


int scotchyyparse (void);



/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_METHODNAME = 3,                 /* METHODNAME  */
  YYSYMBOL_PARAMNAME = 4,                  /* PARAMNAME  */
  YYSYMBOL_VALCASE = 5,                    /* VALCASE  */
  YYSYMBOL_VALDOUBLE = 6,                  /* VALDOUBLE  */
  YYSYMBOL_VALINT = 7,                     /* VALINT  */
  YYSYMBOL_VALSTRING = 8,                  /* VALSTRING  */
  YYSYMBOL_VALSTRAT = 9,                   /* VALSTRAT  */
  YYSYMBOL_VALPARAM = 10,                  /* VALPARAM  */
  YYSYMBOL_VALTEST = 11,                   /* VALTEST  */
  YYSYMBOL_12_ = 12,                       /* '|'  */
  YYSYMBOL_13_ = 13,                       /* '/'  */
  YYSYMBOL_14_ = 14,                       /* '?'  */
  YYSYMBOL_15_ = 15,                       /* ';'  */
  YYSYMBOL_16_ = 16,                       /* ':'  */
  YYSYMBOL_17_ = 17,                       /* '('  */
  YYSYMBOL_18_ = 18,                       /* ')'  */
  YYSYMBOL_19_ = 19,                       /* '{'  */
  YYSYMBOL_20_ = 20,                       /* '}'  */
  YYSYMBOL_21_ = 21,                       /* ','  */
  YYSYMBOL_22_ = 22,                       /* '='  */
  YYSYMBOL_23_ = 23,                       /* '&'  */
  YYSYMBOL_24_ = 24,                       /* '!'  */
  YYSYMBOL_25_ = 25,                       /* '<'  */
  YYSYMBOL_26_ = 26,                       /* '>'  */
  YYSYMBOL_27_ = 27,                       /* '+'  */
  YYSYMBOL_28_ = 28,                       /* '-'  */
  YYSYMBOL_29_ = 29,                       /* '*'  */
  YYSYMBOL_30_ = 30,                       /* '%'  */
  YYSYMBOL_YYACCEPT = 31,                  /* $accept  */
  YYSYMBOL_STRAT = 32,                     /* STRAT  */
  YYSYMBOL_STRATSELECT = 33,               /* STRATSELECT  */
  YYSYMBOL_STRATEMPTY = 34,                /* STRATEMPTY  */
  YYSYMBOL_STRATCONCAT = 35,               /* STRATCONCAT  */
  YYSYMBOL_STRATTEST = 36,                 /* STRATTEST  */
  YYSYMBOL_37_1 = 37,                      /* $@1  */
  YYSYMBOL_38_2 = 38,                      /* $@2  */
  YYSYMBOL_STRATTESTELSE = 39,             /* STRATTESTELSE  */
  YYSYMBOL_STRATGROUP = 40,                /* STRATGROUP  */
  YYSYMBOL_STRATMETHOD = 41,               /* STRATMETHOD  */
  YYSYMBOL_42_3 = 42,                      /* $@3  */
  YYSYMBOL_METHODPARAM = 43,               /* METHODPARAM  */
  YYSYMBOL_44_4 = 44,                      /* $@4  */
  YYSYMBOL_45_5 = 45,                      /* $@5  */
  YYSYMBOL_PARAMLIST = 46,                 /* PARAMLIST  */
  YYSYMBOL_PARAMPARAM = 47,                /* PARAMPARAM  */
  YYSYMBOL_48_6 = 48,                      /* @6  */
  YYSYMBOL_PARAMVAL = 49,                  /* PARAMVAL  */
  YYSYMBOL_50_7 = 50,                      /* @7  */
  YYSYMBOL_TEST = 51,                      /* TEST  */
  YYSYMBOL_TESTOR = 52,                    /* TESTOR  */
  YYSYMBOL_TESTAND = 53,                   /* TESTAND  */
  YYSYMBOL_TESTNOT = 54,                   /* TESTNOT  */
  YYSYMBOL_TESTREL = 55,                   /* TESTREL  */
  YYSYMBOL_TESTRELOP = 56,                 /* TESTRELOP  */
  YYSYMBOL_TESTEXPR1 = 57,                 /* TESTEXPR1  */
  YYSYMBOL_TESTEXPR1OP = 58,               /* TESTEXPR1OP  */
  YYSYMBOL_TESTEXPR2 = 59,                 /* TESTEXPR2  */
  YYSYMBOL_TESTEXPR2OP = 60,               /* TESTEXPR2OP  */
  YYSYMBOL_TESTEXPR3 = 61,                 /* TESTEXPR3  */
  YYSYMBOL_TESTEXPR3OP = 62,               /* TESTEXPR3OP  */
  YYSYMBOL_TESTEXPR4 = 63,                 /* TESTEXPR4  */
  YYSYMBOL_TESTVAL = 64,                   /* TESTVAL  */
  YYSYMBOL_TESTVAR = 65,                   /* TESTVAR  */
  YYSYMBOL_VALSDOUBLE = 66,                /* VALSDOUBLE  */
  YYSYMBOL_VALSINT = 67                    /* VALSINT  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int8 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  13
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   93

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  31
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  37
/* YYNRULES -- Number of rules.  */
#define YYNRULES  65
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  93

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   266


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    24,     2,     2,     2,    30,    23,     2,
      17,    18,    29,    27,    21,    28,     2,    13,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    16,    15,
      25,    22,    26,    14,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    19,    12,    20,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   132,   132,   138,   156,   159,   161,   176,   194,   198,
     202,   198,   226,   229,   234,   239,   243,   247,   246,   308,
     312,   308,   316,   319,   320,   324,   323,   362,   397,   412,
     427,   444,   444,   468,   476,   479,   497,   500,   518,   521,
     537,   541,   544,   563,   567,   571,   577,   593,   596,   600,
     606,   622,   625,   631,   647,   650,   656,   660,   661,   664,
     679,   696,   737,   741,   744,   748
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "METHODNAME",
  "PARAMNAME", "VALCASE", "VALDOUBLE", "VALINT", "VALSTRING", "VALSTRAT",
  "VALPARAM", "VALTEST", "'|'", "'/'", "'?'", "';'", "':'", "'('", "')'",
  "'{'", "'}'", "','", "'='", "'&'", "'!'", "'<'", "'>'", "'+'", "'-'",
  "'*'", "'%'", "$accept", "STRAT", "STRATSELECT", "STRATEMPTY",
  "STRATCONCAT", "STRATTEST", "$@1", "$@2", "STRATTESTELSE", "STRATGROUP",
  "STRATMETHOD", "$@3", "METHODPARAM", "$@4", "$@5", "PARAMLIST",
  "PARAMPARAM", "@6", "PARAMVAL", "@7", "TEST", "TESTOR", "TESTAND",
  "TESTNOT", "TESTREL", "TESTRELOP", "TESTEXPR1", "TESTEXPR1OP",
  "TESTEXPR2", "TESTEXPR2OP", "TESTEXPR3", "TESTEXPR3OP", "TESTEXPR4",
  "TESTVAL", "TESTVAR", "VALSDOUBLE", "VALSINT", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-33)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-32)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int8 yypact[] =
{
       5,   -33,     5,    10,     7,   -33,     5,   -33,     4,   -33,
     -33,    19,    40,   -33,     5,   -33,    27,   -33,    20,   -33,
     -33,   -33,   -33,   -33,    27,    27,   -33,   -33,   -33,    28,
       3,   -33,   -33,    47,    -1,    33,    13,   -33,   -33,   -33,
     -33,   -33,    44,    41,    39,   -33,    42,    27,    27,   -33,
     -33,   -33,    43,    43,   -33,   -33,   -33,    43,   -33,    43,
     -33,    55,   -33,   -33,   -33,     5,     3,   -33,    43,    -4,
      33,    13,   -33,    46,    44,    57,    -9,    14,     8,   -33,
     -33,     5,    63,   -33,   -33,   -33,   -33,     5,   -33,   -33,
       7,   -33,     7
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       6,    17,     6,     0,     2,     4,     5,     8,     0,    12,
      16,    22,     0,     1,     6,     7,     0,    18,     0,    15,
       3,    61,    63,    65,     0,     0,    48,    49,    10,    34,
      36,    38,    41,     0,     0,    47,    51,    54,    57,    58,
      59,    60,     0,     0,     0,    39,     0,     0,     0,    44,
      43,    45,     0,     0,    62,    64,    52,     0,    55,     0,
      25,    20,    24,    40,    56,     6,    35,    37,     0,    42,
      46,    50,    53,     0,     0,     0,    14,     0,     0,    23,
      21,     6,     0,    33,    27,    30,    26,     6,    28,    29,
      13,    11,    32
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -33,   -33,    -2,    66,   -33,    75,   -33,   -33,   -33,   -33,
     -33,   -33,   -33,   -33,   -33,   -33,     9,   -33,   -33,   -33,
     -33,    58,    37,   -21,   -33,   -33,   -22,   -32,    34,   -33,
      29,   -33,    30,   -33,   -33,    12,    15
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     3,     4,     5,     6,     7,     8,    46,    82,     9,
      10,    11,    17,    18,    75,    61,    62,    73,    86,    87,
      28,    29,    30,    31,    32,    52,    33,    34,    35,    57,
      36,    59,    37,    38,    39,    40,    41
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      12,    53,    44,    14,    45,    54,    55,    81,     1,    83,
      13,   -31,    53,    84,    22,    23,    85,    16,    -9,    14,
     -31,   -31,     2,    26,    27,   -31,    48,    67,   -31,   -31,
      69,    21,    64,    22,    23,    26,    27,    53,   -19,    42,
      47,    26,    27,    58,    24,    53,    77,    21,    60,    22,
      23,    25,    14,    47,    26,    27,    65,    64,    19,    63,
      68,    49,    56,    76,    50,    51,    26,    27,    78,    49,
      26,    27,    50,    51,    26,    27,    74,    80,    91,    90,
      20,    15,    43,    79,    66,    92,    71,    70,     0,    72,
      88,     0,     0,    89
};

static const yytype_int8 yycheck[] =
{
       2,    33,    24,    12,    25,     6,     7,    16,     3,     1,
       0,     3,    44,     5,     6,     7,     8,    13,    13,    12,
      12,    13,    17,    27,    28,    17,    23,    48,    20,    21,
      52,     4,    18,     6,     7,    27,    28,    69,    19,    19,
      12,    27,    28,    30,    17,    77,    68,     4,     4,     6,
       7,    24,    12,    12,    27,    28,    14,    18,    18,    18,
      17,    22,    29,    65,    25,    26,    27,    28,    22,    22,
      27,    28,    25,    26,    27,    28,    21,    20,    15,    81,
      14,     6,    24,    74,    47,    87,    57,    53,    -1,    59,
      78,    -1,    -1,    78
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     3,    17,    32,    33,    34,    35,    36,    37,    40,
      41,    42,    33,     0,    12,    36,    13,    43,    44,    18,
      34,     4,     6,     7,    17,    24,    27,    28,    51,    52,
      53,    54,    55,    57,    58,    59,    61,    63,    64,    65,
      66,    67,    19,    52,    57,    54,    38,    12,    23,    22,
      25,    26,    56,    58,     6,     7,    29,    60,    30,    62,
       4,    46,    47,    18,    18,    14,    53,    54,    17,    57,
      59,    61,    63,    48,    21,    45,    33,    57,    22,    47,
      20,    16,    39,     1,     5,     8,    49,    50,    66,    67,
      33,    15,    33
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    31,    32,    33,    33,    34,    34,    35,    35,    37,
      38,    36,    36,    39,    39,    40,    40,    42,    41,    44,
      45,    43,    43,    46,    46,    48,    47,    49,    49,    49,
      49,    50,    49,    49,    51,    52,    52,    53,    53,    54,
      54,    54,    55,    56,    56,    56,    57,    57,    58,    58,
      59,    59,    60,    61,    61,    62,    63,    63,    63,    64,
      64,    65,    66,    66,    67,    67
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     3,     1,     1,     0,     2,     1,     0,
       0,     8,     1,     2,     0,     3,     1,     0,     3,     0,
       0,     5,     0,     3,     1,     0,     4,     1,     1,     1,
       1,     0,     2,     1,     1,     3,     1,     3,     1,     2,
       3,     1,     3,     1,     1,     1,     3,     1,     1,     1,
       3,     1,     1,     3,     1,     1,     3,     1,     1,     1,
       1,     1,     2,     1,     2,     1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      Rf_warning("%s ", Title);                                   \
      Rf_warning(Kind, Value); \
      /*YYFPRINTF (stderr, "\n");*/                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  Rf_warning("Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      Rf_warning(" %d", yybot);
    }
  Rf_warning("\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  Rf_warning("Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      Rf_warning("   $%d = ", yyi + 1);
      Rf_warning(YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      /*YYFPRINTF (stderr, "\n");*/
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  Rf_warning("%s", ("Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  Rf_warning("Entering state %d\n", yystate);
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      Rf_warning("Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      Rf_warning("%s", "Reading a token\n");
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      Rf_warning("%s", "Now at end of input.\n");
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* STRAT: STRATSELECT  */
#line 133 "parser_yy.y"
              {
                parserstratcurr = ((yyvsp[0].STRAT));           /* Save pointer to root of tree */
              }
#line 1369 "y.tab.c"
    break;

  case 3: /* STRATSELECT: STRATSELECT '|' STRATEMPTY  */
#line 139 "parser_yy.y"
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (2)");
                  stratExit  ((yyvsp[-2].STRAT));
                  stratExit  ((yyvsp[0].STRAT));
                  YYABORT;
                }

                strat->tabl                 = parserstrattab;
                strat->type                 = STRATNODESELECT;
                strat->data.select.strat[0] = ((yyvsp[-2].STRAT));
                strat->data.select.strat[1] = ((yyvsp[0].STRAT));

                ((yyval.STRAT)) = strat;
              }
#line 1391 "y.tab.c"
    break;

  case 6: /* STRATEMPTY: %empty  */
#line 161 "parser_yy.y"
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (3)");
                  YYABORT;
                }

                strat->tabl = parserstrattab;
                strat->type = STRATNODEEMPTY;

                ((yyval.STRAT)) = strat;
              }
#line 1409 "y.tab.c"
    break;

  case 7: /* STRATCONCAT: STRATCONCAT STRATTEST  */
#line 177 "parser_yy.y"
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (4)");
                  stratExit  ((yyvsp[-1].STRAT));
                  stratExit  ((yyvsp[0].STRAT));
                  YYABORT;
                }

                strat->tabl                 = parserstrattab;
                strat->type                 = STRATNODECONCAT;
                strat->data.concat.strat[0] = ((yyvsp[-1].STRAT));
                strat->data.concat.strat[1] = ((yyvsp[0].STRAT));

                ((yyval.STRAT)) = strat;
              }
#line 1431 "y.tab.c"
    break;

  case 9: /* $@1: %empty  */
#line 198 "parser_yy.y"
              {
                stratParserSelect (VALTEST);      /* Parse parameter tokens */
              }
#line 1439 "y.tab.c"
    break;

  case 10: /* $@2: %empty  */
#line 202 "parser_yy.y"
              {
                stratParserSelect (VALSTRAT);     /* Parse strategy tokens */
              }
#line 1447 "y.tab.c"
    break;

  case 11: /* STRATTEST: $@1 '/' TEST $@2 '?' STRATSELECT STRATTESTELSE ';'  */
#line 206 "parser_yy.y"
              {
                Strat *           strat;

                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint  ("stratParserParse: out of memory (1)");
                  stratExit ((yyvsp[-2].STRAT));
                  if (((yyvsp[-1].STRAT)) != NULL)
                    stratExit ((yyvsp[-1].STRAT));
                  stratTestExit ((yyvsp[-5].TEST));
                  YYABORT;
                }

                strat->tabl               = parserstrattab;
                strat->type               = STRATNODECOND;
                strat->data.cond.test     = ((yyvsp[-5].TEST));
                strat->data.cond.strat[0] = ((yyvsp[-2].STRAT));
                strat->data.cond.strat[1] = ((yyvsp[-1].STRAT));

                ((yyval.STRAT)) = strat;
              }
#line 1472 "y.tab.c"
    break;

  case 13: /* STRATTESTELSE: ':' STRATSELECT  */
#line 230 "parser_yy.y"
              {
                ((yyval.STRAT)) = ((yyvsp[0].STRAT));
              }
#line 1480 "y.tab.c"
    break;

  case 14: /* STRATTESTELSE: %empty  */
#line 234 "parser_yy.y"
              {
                ((yyval.STRAT)) = NULL;
              }
#line 1488 "y.tab.c"
    break;

  case 15: /* STRATGROUP: '(' STRATSELECT ')'  */
#line 240 "parser_yy.y"
              {
                ((yyval.STRAT)) = ((yyvsp[-1].STRAT));
              }
#line 1496 "y.tab.c"
    break;

  case 17: /* $@3: %empty  */
#line 247 "parser_yy.y"
              {
                Strat *           strat;
                int               meth;
                int               methlen;
                StratMethodTab *  methtab;
                int               i, j;

                meth    =
                methlen = 0;                      /* No method recognized yet   */
                methtab = parserstrattab->methtab; /* Point to the method table */
                for (i = 0; methtab[i].name != NULL; i ++) {
                  if ((strncasecmp (((yyvsp[0].STRING)),         /* Find longest matching code name */
                       methtab[i].name,
                       j = strlen (methtab[i].name)) == 0) &&
                      (j > methlen)) {
                    meth    = methtab[i].meth;
                    methlen = j;
                  }
                }
                if (methlen == 0) {               /* If method name not known */
                  errorPrint ("stratParserParse: invalid method name \"%s\", before \"%s\"",
                              ((yyvsp[0].STRING)), stratParserRemain ());
                  YYABORT;
                }
                if ((strat = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (5)");
                  YYABORT;
                }

                strat->tabl             = parserstrattab;
                strat->type             = STRATNODEMETHOD;
                strat->data.method.meth = meth;   /* Set method type         */
                if (methtab[meth].data != NULL)   /* If default values exist */
                  memcpy (&strat->data.method.data, /* Set values to default */
                          methtab[meth].data,
                          sizeof (StratNodeMethodData));

                parserstratcurr = strat;          /* Structure available for parameter processing */
              }
#line 1540 "y.tab.c"
    break;

  case 18: /* STRATMETHOD: METHODNAME $@3 METHODPARAM  */
#line 287 "parser_yy.y"
              {
                StratParamTab *   paratab;
                int               i;

                paratab = parserstrattab->paratab; /* Point to the parameter table */
                for (i = 0; paratab[i].name != NULL; i ++) {
                  if ((paratab[i].meth == parserstratcurr->data.method.meth) && /* If a strategy parameter found for this method */
                      (paratab[i].type == STRATPARAMSTRAT)) {
                    if (*((Strat **) ((byte *) &parserstratcurr->data.method.data + /* And this parameter has not been set */
                        (paratab[i].dataofft - paratab[i].database))) == NULL)
                      errorPrintW ("stratParserParse: strategy parameter \"%s\" of method \"%s\" not set, before \"%s\"",
                                   paratab[i].name, parserstrattab->methtab[parserstratcurr->data.method.meth].name, stratParserRemain ());
                  }
                }

                ((yyval.STRAT)) = parserstratcurr;           /* Return current structure */
                parserstratcurr = NULL;           /* No current structure     */
              }
#line 1563 "y.tab.c"
    break;

  case 19: /* $@4: %empty  */
#line 308 "parser_yy.y"
              {
                stratParserSelect (VALPARAM);     /* Parse parameter tokens */
              }
#line 1571 "y.tab.c"
    break;

  case 20: /* $@5: %empty  */
#line 312 "parser_yy.y"
              {
                stratParserSelect (VALSTRAT);     /* Parse strategy tokens */
              }
#line 1579 "y.tab.c"
    break;

  case 25: /* @6: %empty  */
#line 324 "parser_yy.y"
              {
                int               para;
                int               paralen;
                StratParamTab *   paratab;
                int               i, j;

                para    =
                paralen = 0;                      /* No parameter recognized yet   */
                paratab = parserstrattab->paratab; /* Point to the parameter table */
                for (i = 0; paratab[i].name != NULL; i ++) {
                  if ((paratab[i].meth == parserstratcurr->data.method.meth) &&
                      (strncasecmp (((yyvsp[0].STRING)),         /* Find longest matching parameter name */
                                    paratab[i].name,
                                    j = strlen (paratab[i].name)) == 0) &&
                      (j > paralen)) {
                    para    = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid method parameter name \"%s\", before \"%s\"",
                              ((yyvsp[0].STRING)), stratParserRemain ());
                  YYABORT;
                }

                ((yyval.SAVE)).tabl = parserstrattab; /* Save current strategy tables */
                parserparamcurr = &paratab[para]; /* Save current parameter value */
                stratParserSelect (parsermethtokentab[parserparamcurr->type & ~STRATPARAMDEPRECATED]); /* Get non-deprecated type */
                if (parserparamcurr->type == STRATPARAMSTRAT) /* If parameter is a strategy           */
                  parserstrattab = (StratTab *) parserparamcurr->datasltr; /* Use new strategy tables */
              }
#line 1615 "y.tab.c"
    break;

  case 26: /* PARAMPARAM: PARAMNAME @6 '=' PARAMVAL  */
#line 356 "parser_yy.y"
              {
                stratParserSelect (VALPARAM);     /* Go-on reading parameters        */
                parserstrattab = ((yyvsp[-2].SAVE)).tabl; /* Restore current strategy tables */
              }
#line 1624 "y.tab.c"
    break;

  case 27: /* PARAMVAL: VALCASE  */
#line 363 "parser_yy.y"
              {
                char              c;              /* Character read             */
                char *            p;              /* Pointer to selector string */
                int               i;              /* Index in selector string   */

                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
                  c = ((yyvsp[0].CASEVAL));                       /* First, use char as is                                   */
                  for (p = (char *) parserparamcurr->datasltr, i = 0;
                       (*p != '\0') && (*p != c);
                       p ++, i ++) ;
                  if (*p == '\0') {               /* Char was not found         */
                    c = tolower (c);              /* Convert char to lower case */
                    for (p = (char *) parserparamcurr->datasltr, i = 0;
                         (*p != '\0') && (*p != c);
                         p ++, i ++) ;
                    if (*p == '\0') {
                      errorPrint ("stratParserParse: invalid method parameter switch \"%s=%c\", before \"%s\"",
                                  parserparamcurr->name, ((yyvsp[0].CASEVAL)), stratParserRemain ());
                      YYABORT;
                    }
                  }

#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (int)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (1)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((int *) ((byte *) &parserstratcurr->data.method.data +
                             (parserparamcurr->dataofft -
                              parserparamcurr->database))) = i;
                }
              }
#line 1663 "y.tab.c"
    break;

  case 28: /* PARAMVAL: VALSDOUBLE  */
#line 398 "parser_yy.y"
              {
                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (double)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (2)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((double *) ((byte *) &parserstratcurr->data.method.data +
                                (parserparamcurr->dataofft -
                                 parserparamcurr->database))) = ((yyvsp[0].DOUBLE));
                }
              }
#line 1682 "y.tab.c"
    break;

  case 29: /* PARAMVAL: VALSINT  */
#line 413 "parser_yy.y"
              {
                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (INT)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (3)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((INT *) ((byte *) &parserstratcurr->data.method.data +
                             (parserparamcurr->dataofft -
                              parserparamcurr->database))) = (INT) ((yyvsp[0].INTEGER));
                }
              }
#line 1701 "y.tab.c"
    break;

  case 30: /* PARAMVAL: VALSTRING  */
#line 428 "parser_yy.y"
              {
                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + strlen ((yyvsp[0].STRING)) + 1) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (4)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  strcpy ((char *) ((byte *) &parserstratcurr->data.method.data +
                                    (parserparamcurr->dataofft -
                                     parserparamcurr->database)),
                          ((yyvsp[0].STRING)));
                }
              }
#line 1721 "y.tab.c"
    break;

  case 31: /* @7: %empty  */
#line 444 "parser_yy.y"
              {
                ((yyval.SAVE)).strat = parserstratcurr;
                ((yyval.SAVE)).param = parserparamcurr;
                parserstratcurr  = NULL;
                parserparamcurr  = NULL;
              }
#line 1732 "y.tab.c"
    break;

  case 32: /* PARAMVAL: @7 STRATSELECT  */
#line 451 "parser_yy.y"
              {
                parserstratcurr = ((yyvsp[-1].SAVE)).strat; /* Restore current method    */
                parserparamcurr = ((yyvsp[-1].SAVE)).param; /* Restore current parameter */

                if ((parserparamcurr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((parserparamcurr->dataofft - parserparamcurr->database + sizeof (Strat *)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (5)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((Strat **) ((byte *) &parserstratcurr->data.method.data +
                                (parserparamcurr->dataofft -
                                 parserparamcurr->database))) = ((yyvsp[0].STRAT));
                }
              }
#line 1754 "y.tab.c"
    break;

  case 33: /* PARAMVAL: error  */
#line 469 "parser_yy.y"
              {
                errorPrint ("stratParserParse: invalid value for parameter \"%s\" of method \"%s\", before \"%s\"",
                            parserparamcurr->name, parserstratcurr->tabl->methtab[parserstratcurr->data.method.meth].name, stratParserRemain ());
                YYABORT;
              }
#line 1764 "y.tab.c"
    break;

  case 35: /* TESTOR: TESTOR '|' TESTAND  */
#line 480 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (6)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }

                test->typetest     = STRATTESTOR;
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ((yyvsp[-2].TEST));
                test->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = test;
              }
#line 1786 "y.tab.c"
    break;

  case 37: /* TESTAND: TESTAND '&' TESTNOT  */
#line 501 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (7)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }

                test->typetest     = STRATTESTAND;
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ((yyvsp[-2].TEST));
                test->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = test;
              }
#line 1808 "y.tab.c"
    break;

  case 39: /* TESTNOT: '!' TESTNOT  */
#line 522 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (8)");
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }

                test->typetest     = STRATTESTNOT;
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = test;
              }
#line 1828 "y.tab.c"
    break;

  case 40: /* TESTNOT: '(' TESTOR ')'  */
#line 538 "parser_yy.y"
              {
                ((yyval.TEST)) = ((yyvsp[-1].TEST));
              }
#line 1836 "y.tab.c"
    break;

  case 42: /* TESTREL: TESTEXPR1 TESTRELOP TESTEXPR1  */
#line 545 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (9)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }
                test->typetest     = ((yyvsp[-1].TESTOP));
                test->typenode     = STRATPARAMLOG;
                test->data.test[0] = ((yyvsp[-2].TEST));
                test->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = test;
              }
#line 1857 "y.tab.c"
    break;

  case 43: /* TESTRELOP: '<'  */
#line 564 "parser_yy.y"
              {
                ((yyval.TESTOP)) = STRATTESTLT;
              }
#line 1865 "y.tab.c"
    break;

  case 44: /* TESTRELOP: '='  */
#line 568 "parser_yy.y"
              {
                ((yyval.TESTOP)) = STRATTESTEQ;
              }
#line 1873 "y.tab.c"
    break;

  case 45: /* TESTRELOP: '>'  */
#line 572 "parser_yy.y"
              {
                ((yyval.TESTOP)) = STRATTESTGT;
              }
#line 1881 "y.tab.c"
    break;

  case 46: /* TESTEXPR1: TESTEXPR1 TESTEXPR1OP TESTEXPR2  */
#line 578 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (10)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }
                test->typetest     = ((yyvsp[-1].TESTOP));
                test->data.test[0] = ((yyvsp[-2].TEST));
                test->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = test;
              }
#line 1901 "y.tab.c"
    break;

  case 48: /* TESTEXPR1OP: '+'  */
#line 597 "parser_yy.y"
              {
                ((yyval.TESTOP)) = STRATTESTADD;
              }
#line 1909 "y.tab.c"
    break;

  case 49: /* TESTEXPR1OP: '-'  */
#line 601 "parser_yy.y"
              {
                ((yyval.TESTOP)) = STRATTESTSUB;
              }
#line 1917 "y.tab.c"
    break;

  case 50: /* TESTEXPR2: TESTEXPR2 TESTEXPR2OP TESTEXPR3  */
#line 607 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  errorPrint    ("stratParserParse: out of memory (11)");
                  YYABORT;
                }
                test->typetest     = ((yyvsp[-1].TESTOP));
                test->data.test[0] = ((yyvsp[-2].TEST));
                test->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = test;
              }
#line 1937 "y.tab.c"
    break;

  case 52: /* TESTEXPR2OP: '*'  */
#line 626 "parser_yy.y"
              {
                ((yyval.TESTOP)) = STRATTESTMUL;
              }
#line 1945 "y.tab.c"
    break;

  case 53: /* TESTEXPR3: TESTEXPR3 TESTEXPR3OP TESTEXPR4  */
#line 632 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (12)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }
                test->typetest     = ((yyvsp[-1].TESTOP));
                test->data.test[0] = ((yyvsp[-2].TEST));
                test->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = test;
              }
#line 1965 "y.tab.c"
    break;

  case 55: /* TESTEXPR3OP: '%'  */
#line 651 "parser_yy.y"
              {
                ((yyval.TESTOP)) = STRATTESTMOD;
              }
#line 1973 "y.tab.c"
    break;

  case 56: /* TESTEXPR4: '(' TESTEXPR1 ')'  */
#line 657 "parser_yy.y"
              {
                ((yyval.TEST)) = ((yyvsp[-1].TEST));
              }
#line 1981 "y.tab.c"
    break;

  case 59: /* TESTVAL: VALSDOUBLE  */
#line 665 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (13)");
                  YYABORT;
                }

                test->typetest        = STRATTESTVAL;
                test->typenode        = STRATPARAMDOUBLE;
                test->data.val.valdbl = ((yyvsp[0].DOUBLE));

                ((yyval.TEST)) = test;
              }
#line 2000 "y.tab.c"
    break;

  case 60: /* TESTVAL: VALSINT  */
#line 680 "parser_yy.y"
              {
                StratTest *       test;

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (14)");
                  YYABORT;
                }

                test->typetest        = STRATTESTVAL;
                test->typenode        = STRATPARAMINT;
                test->data.val.valint = ((yyvsp[0].INTEGER));

                ((yyval.TEST)) = test;
              }
#line 2019 "y.tab.c"
    break;

  case 61: /* TESTVAR: PARAMNAME  */
#line 697 "parser_yy.y"
              {
                StratTest *       test;
                StratParamTab *   condtab;
                int               para;
                int               paralen;
                int               i, j;

                para    =
                paralen = 0;                      /* No parameter recognized yet */
                condtab = parserstrattab->condtab; /* Point to parameter table   */
                for (i = 0; condtab[i].name != NULL; i ++) {
                  if ((strncasecmp (((yyvsp[0].STRING)),         /* Find longest matching parameter name */
                                    condtab[i].name,
                                    j = strlen (condtab[i].name)) == 0) &&
                      (j > paralen)) {
                    para    = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid graph parameter name \"%s\", before \"%s\"",
                              ((yyvsp[0].STRING)), stratParserRemain ());
                  YYABORT;
                }

                if ((test = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (15)");
                  YYABORT;
                }

                test->typetest          = STRATTESTVAR;
                test->typenode          = condtab[para].type;
                test->data.var.datatab  = parserstrattab;
                test->data.var.datadisp = condtab[para].dataofft -
                                          condtab[para].database;

                ((yyval.TEST)) = test;
              }
#line 2062 "y.tab.c"
    break;

  case 62: /* VALSDOUBLE: TESTEXPR1OP VALDOUBLE  */
#line 738 "parser_yy.y"
              {
                ((yyval.DOUBLE)) = (((yyvsp[-1].TESTOP)) == STRATTESTSUB) ? - ((yyvsp[0].DOUBLE)) : ((yyvsp[0].DOUBLE));
              }
#line 2070 "y.tab.c"
    break;

  case 64: /* VALSINT: TESTEXPR1OP VALINT  */
#line 745 "parser_yy.y"
              {
                ((yyval.INTEGER)) = (((yyvsp[-1].TESTOP)) == STRATTESTSUB) ? - ((yyvsp[0].INTEGER)) : ((yyvsp[0].INTEGER));
              }
#line 2078 "y.tab.c"
    break;


#line 2082 "y.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 751 "parser_yy.y"


/*
**  The static and global definitions (bis).
**  These are put at the end of the file because
**  the token values that they use are not yet
**  defined in the first section of the file.
*/

unsigned int                parsermethtokentab[] = { /* Table for parameter/token type conversion */
                              VALCASE,
                              VALDOUBLE,
                              VALINT,
                              -1,                 /* No logical parameters */
                              VALSTRAT,
                              VALSTRING,
                              -1                  /* One more value to detect array overflow */
                            };

/************************************/
/*                                  */
/* These routines drive the parser. */
/*                                  */
/************************************/

/* This routine is the entry point for
** the strategy parser.
** It returns:
** - !NULL  : pointer to the strategy.
** - NULL   : on error.
*/

Strat *
stratParserParse (
const StratTab * const      strattab,             /*+ Pointer to parsing tables +*/
const char * const          string)               /*+ Strategy string to parse  +*/
{
  yyclearin;                                      /* Reset the parser state */

#ifdef SCOTCH_DEBUG_PARSER3
  yydebug = 1;                                    /* Set debugging if needed */
#endif /* SCOTCH_DEBUG_PARSER3 */

  stratParserInit (string);                       /* Initialize the lexical parser           */
  parserstrattab  = strattab;                     /* Point to the parsing tables             */
  parserstratcurr = NULL;                         /* Clear up the temporary strategy pointer */

  if (yyparse () != 0) {                          /* Parse the strategy string */
    if (parserstratcurr != NULL)
      stratExit (parserstratcurr);
    return (NULL);
  }

  return (parserstratcurr);                       /* Return strategy pointer */
}

/* This routine displays the parser error message.
** It returns:
** - 1  : in all cases.
*/

static
int
yyerror (
const char * const          errstr)
{
  errorPrint ("stratParserParse: invalid strategy string, before \"%s\"", stratParserRemain ());
  return     (1);
}
