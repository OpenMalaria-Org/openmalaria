/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef OM_util_errors
#define OM_util_errors

#include <stdexcept>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

/* Macros to ease use of traced_exception.
 * msg should be obvious, and code is the program exit code.
 */
#define TRACED_EXCEPTION( msg, code ) ::OM::util::traced_exception( (msg), __FILE__, __LINE__, (code) )
#define TRACED_EXCEPTION_DEFAULT( msg ) ::OM::util::traced_exception( (msg), __FILE__, __LINE__ )
#define SWITCH_DEFAULT_EXCEPTION ::OM::util::traced_exception( OM::util::Messages::SwitchDefault, __FILE__, __LINE__, ::OM::util::Error::SwitchDefault )

// #define OM_NO_STACK_TRACE

#ifdef OM_NO_STACK_TRACE
// do nothing
#elif defined __GNU_LIBRARY__
#define OM_GNU_STACK_TRACE
#endif

// Macro to make checking XML inputs easier
#define XML_ASSERT( cond, msg ) if(!(cond)) throw ::OM::util::xml_scenario_error(msg)

#define OM_ASSERT_EQ( a, b ) ::OM::util::do_assertion( (a), (b), __FILE__, __LINE__ )
#define OM_ASSERT_EQ_TOL( a, b, tol ) ::OM::util::do_assertion_tol( (a), (b), (tol), __FILE__, __LINE__ )

/** Standard exception classes for OpenMalaria. */
namespace OM { namespace util {
    
    /** Different exit codes which may be used.
     * 
     * These are intended to categorise errors, not diagnose them. Always use
     * error messages or stack traces to diagnose errors. (The only reason we
     * don't just return 1 is that BOINC uses the codes to correlate errors.)
     * 
     * We start exit codes at 64, as in /usr/include/sysexits.h; as far as I am
     * aware, linux is the only platform targetted which has such narrow ranges
     * of exit codes available (hence these codes should work on Windows and
     * Mac OS X too). */
    namespace Error { enum ErrorCodes {
        /// no error (exit code 0)
        None = 0,
        Default = 64,
        TracedDefault,
        Assertion,
        /// any error from Code Synthesis's XSD
        XSD,
        /// any checkpointing error (these don't occur often enough to need segregating here)
        Checkpoint,
        /// invalid scenario file
        XmlScenario,
        GSL,
        FileExists,		// wanted to create a file but it already exists!
        FileIO,			// any other file read/write error
        EffectiveEIR,
        NumNewInfections,
        InitialKappa,
        VectorWarmup,
        CommandLine,
        SumWeight,
        VectorFitting,
        InfLambda,
        NotImplemented,
        WHFeatures,     // these are probably all code errors
        SwitchDefault,
        InputResource,
        PkPd,
        NoStartDate,
        Max
    }; }
    namespace Messages {
        extern const char *SwitchDefault;
    }
    
    // As in the "Advanced Bash-Scripting Guide"; not directly relevant to C++
    // but gives some idea what codes make sense to use.
    static_assert( Error::Max <= 113, "Error::Max <= 113" );
    
    /** Base OpenMalaria exception class.
     * 
     * Throw this directly for an error not fitting one of the below
     * subclasses.
     * 
     * When throwing this or a subclass, give a message which should be
     * informative and clear from a user-perspective, unless the result is a
     * code error. In this case use traced_exception. */
    class base_exception : public runtime_error {
    public:
        // for common errors, specify a unique TRACED_DEFAULT to help segregate
        // errors in BOINC server summaries; otherwise leave the default
        explicit base_exception(const string& msg, int code=Error::Default);
        /** Get the classification. This should be some value of the
         * Error::ErrorCodes enum. */
        inline int getCode() const{
            return errCode;
        }
        /** Get the message. Use this as opposed to what() to enable overriding
         * (without dealing with the messy details of how exactly each standard
         * library defines what()). */
        virtual const char* message() const{
            return what();
        }
    private:
        int errCode;
    };
    
    /** Extension of base_exception for not-yet-implemented features. */
    class unimplemented_exception : public base_exception {
    public:
        explicit unimplemented_exception(const string& msg);
        
        /** Complete the error message.
         *
         * Warning: this function is not very robust. References to the
         * returned result should not survive beyond a second call to this
         * function (from any unimplemented_exception object). */
        virtual const char* message() const;
    };
    
    /** Extension of base_exception which tries to get a stack trace. Also
     * prints a message indicating that this is likely a code error.
     * 
     * Should be used when catching code errors to get extra debugging
     * information. */
    class traced_exception : public base_exception {
    public:
        /** Create a stack-trace and store in this object.
         * 
         * @param msg Error message
         * @param file File in which it occurred (may pass nullptr)
         * @param line Line of error (may pass 0)
         * @param start Index of first stack frame of interest. 0 is this
         * constructor, 1 the code creating this exception (usually the location
         * of the problem, but sometimes an error handler), 2 the next frame etc.
         * @param code Exit code.
         * 
         * See TRACED_EXCEPTION and TRACED_EXCEPTION_DEFAULT macros which
         * reduce the number of arguments which need to be passed.
         */
        explicit traced_exception(const string& msg, const char *file, int line, int code=Error::TracedDefault, int start=1);
        virtual ~traced_exception() throw();
    private:
        const char *file;
        int line, start;
#ifdef OM_GNU_STACK_TRACE
        size_t length;
        char** trace;
#endif
        friend ostream& operator<<(ostream& stream, const traced_exception& e);
    };
    /// Print trace to stream (a series of lines, each with \\n appended)
    ostream& operator<<(ostream& stream, const traced_exception& e);
    
    /** Thrown to indicate an error in the scenario.xml file.
     * 
     * Since such errors are always user-input errors, collecting a stack
     * trace doesn't help. Just provide a useful message! */
    class xml_scenario_error : public base_exception {
    public:
        explicit xml_scenario_error(const string& msg);
    };
    
    /** Thrown to indicate a badly formatted parameter. */
    class format_error : public xml_scenario_error {
    public:
        explicit format_error (const string& what_arg);
    };
    
    /** Thrown to indicate an error while loading/saving a checkpoint.
    *
    * Prepends "Error reading checkpoint: " to the message. */
    class checkpoint_error : public traced_exception {
    public:
        explicit checkpoint_error(const string& msg);
    };
    
    /** Thrown due to command-line error or to halt, when a command-line
     * argument prompts an early exit. Not an error when errCode==0. */
    class cmd_exception : public base_exception {
    public:
        explicit cmd_exception(const string& msg, int code=Error::CommandLine);
    };
    
    /** Set up a GSL handler to throw a traced exception. */
    void set_gsl_handler();
    
    template<typename A, typename B>
    void do_assertion(A a, B b, const char *file, int line) {
        if( a != b ) {
                ostringstream oss;
                oss << "Expected " << a << " == " << b;
                throw traced_exception( oss.str(), file, line, Error::Assertion );
        }
    }
    
    template<typename T>
    void do_assertion_tol(T a, T b, T tol, const char *file, int line) {
        if( !((a < b + tol ) && (b < a + tol)) ){
                ostringstream oss;
                oss << "Expected " << a << " == " << b << " to tolerance " << tol;
                throw traced_exception( oss.str(), file, line, Error::Assertion );
        }
    }
} }
#endif
