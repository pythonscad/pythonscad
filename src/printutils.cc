#include "printutils.h"
#include <sstream>
#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/filesystem.hpp>
#include "exceptions.h"

namespace fs = boost::filesystem;

std::list<std::string> print_messages_stack;
std::list<struct Message> log_messages_stack;
OutputHandlerFunc *outputhandler = nullptr;
void *outputhandler_data = nullptr;
OutputHandlerFunc2 *outputhandler2 = nullptr;
void *outputhandler_data2 = nullptr;
std::string OpenSCAD::debug("");
bool OpenSCAD::quiet = false;
bool OpenSCAD::hardwarnings = false;
bool OpenSCAD::parameterCheck = true;
bool OpenSCAD::rangeCheck = false;

boost::circular_buffer<std::string> lastmessages(5);
boost::circular_buffer<struct Message> lastlogmessages(5);

int count=0;

namespace {
	bool no_throw;
	bool deferred;
}

void set_output_handler(OutputHandlerFunc *newhandler, void *userdata)
{
	outputhandler = newhandler;
	outputhandler_data = userdata;
}

void set_output_handler2(OutputHandlerFunc2 *newhandler, void *userdata)
{
	outputhandler2 = newhandler;
	outputhandler_data2 = userdata;
}

void no_exceptions_for_warnings()
{
	no_throw = true;
	deferred = false;
}

bool would_have_thrown()
{
    const auto would_throw = deferred;
    no_throw = false;
    deferred = false;
    return would_throw;
}

void print_messages_push()
{
	print_messages_stack.push_back(std::string());
}

void print_messages_pop()
{
	std::string msg = print_messages_stack.back();
	print_messages_stack.pop_back();
	if (print_messages_stack.size() > 0 && !msg.empty()) {
		if (!print_messages_stack.back().empty()) {
			print_messages_stack.back() += "\n";
		}
		print_messages_stack.back() += msg;
	}
}

void PRINT(const std::string &msg)
{
	if (msg.empty()) return;
	if (print_messages_stack.size() > 0) {
		if (!print_messages_stack.back().empty()) {
			print_messages_stack.back() += "\n";
		}
		print_messages_stack.back() += msg;
	}
	PRINT_NOCACHE(msg);
}

void PRINT_NOCACHE(const std::string &msg)
{
	if (msg.empty()) return;

	if (boost::starts_with(msg, "WARNING") || boost::starts_with(msg, "ERROR") || boost::starts_with(msg, "TRACE")) {
		size_t i;
		for (i=0;i<lastmessages.size();i++) {
			if (lastmessages[i] != msg) break;
		}
		if (i == 5) return; // Suppress output after 5 equal ERROR or WARNING outputs.
		else lastmessages.push_back(msg);
	}
	if(!deferred)
		if (!OpenSCAD::quiet || boost::starts_with(msg, "ERROR")) {
			if (!outputhandler) {
				fprintf(stderr, "%s\n", msg.c_str());
			} else {
				outputhandler(msg, outputhandler_data);
			}
		}
	if(!std::current_exception()) {
		if((OpenSCAD::hardwarnings && boost::starts_with(msg, "WARNING")) || (no_throw && boost::starts_with(msg, "ERROR"))){
			if(no_throw)
				deferred = true;
			else
				throw HardWarningException(msg);
		}
	}
}

void LOG(const Message &msg)
{
	if (!outputhandler2) {
		// fprintf(stderr, "%s\n", msg.c_str());
		} else {
			outputhandler2(msg, outputhandler_data2);
		}
}

void PRINTDEBUG(const std::string &filename, const std::string &msg)
{
	// see printutils.h for usage instructions
	if (OpenSCAD::debug=="") return;
	std::string shortfname = fs::path(filename).stem().generic_string();
	std::string lowshortfname(shortfname);
	boost::algorithm::to_lower(lowshortfname);
	std::string lowdebug(OpenSCAD::debug);
	boost::algorithm::to_lower(lowdebug);
	if (OpenSCAD::debug=="all" ||
			lowdebug.find(lowshortfname) != std::string::npos) {
		PRINT_NOCACHE( shortfname+": "+ msg );
	}
}

const std::string& quoted_string(const std::string& str)
{
	static std::string buf;
	buf = str;
	boost::replace_all(buf, "\n", "\\n");
	return buf;
}

std::string two_digit_exp_format( std::string doublestr )
{
#ifdef _WIN32
	size_t exppos = doublestr.find('e');
	if ( exppos != std::string::npos) {
		exppos += 2;
		if ( doublestr[exppos] == '0' ) doublestr.erase(exppos,1);
	}
#endif
	return doublestr;
}

std::string two_digit_exp_format(double x)
{
	return two_digit_exp_format(std::to_string(x));
}

#include <set>

std::set<std::string> printedDeprecations;

void printDeprecation(const std::string &str)
{
	if (printedDeprecations.find(str) == printedDeprecations.end()) {
		printedDeprecations.insert(str);
		std::string msg = "DEPRECATED: " + str;
		PRINT(msg);
	}
}

void resetSuppressedMessages()
{
	printedDeprecations.clear();
	lastmessages.clear();
}
