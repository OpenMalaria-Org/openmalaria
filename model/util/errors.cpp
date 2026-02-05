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

#include "util/errors.h"
#include <ostream>
#include <iostream>
#include <gsl_errno.h>
#include <cstring>

#ifdef __GNU_LIBRARY__
#include <execinfo.h>
#include <stdlib.h>
#include <cxxabi.h>
#endif

/** Standard exception classes for OpenMalaria. */
namespace OM { namespace util {

const char *Messages::SwitchDefault = "no valid option selected in a switch statement or if/else tree";

base_exception::base_exception(const string& msg, int code) :
    runtime_error(msg),
    errCode(code)
{}

traced_exception::traced_exception(const string& msg, const char *f, int l, int code, int s) :
    base_exception(msg,code),
    file(f), line(l), start(s)
{
#ifdef OM_GNU_STACK_TRACE
    // http://www.gnu.org/software/libc/manual/html_node/Backtraces.html
// probably 100 is big enough
#define MAX_STACK_SIZE 100
    void *array[MAX_STACK_SIZE];
    
    length = backtrace (array, MAX_STACK_SIZE);
    trace = backtrace_symbols (array, length);
#endif
}
traced_exception::~traced_exception() throw(){
#ifdef OM_GNU_STACK_TRACE
    free (trace);
#endif
}
ostream& operator<<(ostream& stream, const traced_exception& e){
    stream<<"Call stack";
    if( e.file != 0 )
        stream << ", starting from "<<e.file<<':'<<e.line;
    stream<<":\n";
#ifdef OM_NO_STACK_TRACE
    stream << "backtrace capture disabled\n";
#elif defined OM_GNU_STACK_TRACE
    // demangle output: http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
    // spec: http://www.ib.cnea.gov.ar/~oop/biblio/libstdc++/namespaceabi.html
    
    //string binary;
    // Usually e.start=1: we skip the first frame (traced_exception constructor)
    for( size_t i=e.start; i<e.length; ++i ){
        char *line = e.trace[i];
        char *lb = strchr(line, '(');
        char *plus = 0;
        char *rb = 0;
        char *demangled = 0;
        if( lb != 0 ){
            /* Print out the path to the binary? Seems pointless.
            if( binary != string( line, lb-line ) ){
                binary = string( line, lb-line );
                stream << "in " << binary << ":\n";
            }
            */
            *lb = '\0';
            plus = strchr(lb+1,'+');
            if( plus != 0 ){
                *plus = '\0';
                int status;
                demangled = abi::__cxa_demangle(lb+1, 0, 0, &status);
                /*if( status != 0 ){
                    cerr << "status: "<<status<<"\tstring:"<<lb+1<<endl;
                }*/
                rb = strchr(plus+1, ')');
            }else{
                rb = strchr(lb+1, ')');
            }
            if( rb != 0 ){
                *rb = '\0';
            }
        }
        if( lb == 0 ){
            stream << line;
        }else{
            if( plus != 0 ){
                stream << "+" << (plus+1);
            }
            if( demangled != 0 ){
                stream << '\t' << demangled;
                free(demangled);
            }else if( lb != 0 && rb > lb ){
                stream << '\t' << (lb+1) << "()";
            }
        }
        stream << '\n';
    }
#else
    stream << "sorry, no trace from this platform!\n";
#endif
    return stream;
}

unimplemented_exception::unimplemented_exception(const string& msg) :
        base_exception(msg, Error::NotImplemented)
{}

const string msg_pre = "Sorry, this feature hasn't been implemented by OpenMalaria yet: ";
vector<char> msg_buf;

const char *unimplemented_exception::message() const{
    const char *base_msg = base_exception::what();
    size_t base_len = strlen( base_msg );
    size_t required = base_len + msg_pre.size();
    if( msg_buf.size() == 0 ){
        msg_buf.assign( required, '0' );
        memcpy( &msg_buf[0], msg_pre.c_str(), msg_pre.size() );
    }else if( msg_buf.size() < required ){
        msg_buf.resize( required );
    }
    memcpy( &msg_buf[msg_pre.size()], base_msg, base_len );
    return &msg_buf[0];
}

xml_scenario_error::xml_scenario_error(const string& msg) :
    base_exception(msg, Error::XmlScenario)
{}

format_error::format_error(const string& msg) :
    xml_scenario_error(msg)
{}

checkpoint_error::checkpoint_error(const string& msg) :
    traced_exception(msg, 0, 0, Error::Checkpoint)
{}

cmd_exception::cmd_exception(const string& msg, int code) :
    base_exception(msg, code)
{}

void gsl_handler( const char *reason,
                  const char *file,
                  int line,
                  int gsl_errno){
    throw traced_exception(reason, file, line, Error::GSL, 2);
}

void set_gsl_handler() {
    gsl_set_error_handler( &gsl_handler );
}
   
} }
