#include "genlang.h"

std::vector<std::string> mapping_name;
std::vector<std::string> mapping_code;
std::vector<int> mapping_level;

std::vector<std::shared_ptr<AbstractNode>> shows;
std::shared_ptr<AbstractNode> genlang_result_node = nullptr;
int language=LANG_SCAD;

LanguageDesc lang_scad("scad",
		[](const char *fname, const char *content){
		  return boost::algorithm::ends_with(fname, ".scad")?1:0;
		  }
		);
// LanguageDesc lang_js({.suffix="js"});

extern LanguageDesc lang_scad;
extern LanguageDesc lang_py;
extern LanguageDesc lang_js;
extern LanguageDesc lang_lua;

LanguageDesc *languageDesc[4]=
{
	&lang_scad,
#ifdef ENABLE_PYTHON	
	&lang_py,
#else
	nullptr,
#endif	
		
#ifdef ENABLE_JS	
	&lang_js,
#else
	nullptr,
#endif	
		
#ifdef ENABLE_LUA	
	&lang_lua
#else
	nullptr
#endif	
		
};

void show_final(void)
{
  mapping_name.clear(); 
  mapping_code.clear();
  mapping_level.clear();
  if(shows.size() == 1) genlang_result_node = shows[0];
  else {
    DECLARE_INSTANCE
    genlang_result_node = std::make_shared<CsgOpNode>(instance, OpenSCADOperator::UNION); 
    genlang_result_node -> children = shows;
  }
  shows.clear();

}


