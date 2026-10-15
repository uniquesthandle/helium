CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++20 -Wall -Wextra
CORE      = he_api.cpp he_eos.cpp he_table.cpp

all: he_table_build he_calc he_test he_fit

he_table_build: he_table_build.cpp he_ref.cpp $(CORE) he_api.h he_eos.h he_table.h he_ref.h he_console.h
	$(CXX) $(CXXFLAGS) he_table_build.cpp he_ref.cpp $(CORE) -o $@

#Построитель с эталоном CoolProp (пути укажите свои):
#   make coolprop CP_INC="/e/misc/my_coolprop_project/externals/CoolProp/include" CP_LIB="/e/misc/my_coolprop_project/build/externals/CoolProp"

CP_INC = /c/code/my_coolprop_project/externals/CoolProp/include -I/c/code/my_coolprop_project/build/_deps/fmt-src/include
CP_LIB = /c/code/my_coolprop_project/build/CoolProp/Release

coolprop: he_table_build.cpp he_ref.cpp $(CORE)
	$(CXX) $(CXXFLAGS) -DUSE_COOLPROP -I$(CP_INC) he_table_build.cpp he_ref.cpp $(CORE) -L$(CP_LIB) -lCoolProp -o he_table_build

he_calc: he_calc.cpp $(CORE) he_api.h he_eos.h he_table.h
	$(CXX) $(CXXFLAGS) he_calc.cpp $(CORE) -o $@

he_test: he_test.cpp $(CORE) he_api.h he_eos.h he_table.h
	$(CXX) $(CXXFLAGS) he_test.cpp $(CORE) -o $@

he_fit: he_fit.cpp he_eos.cpp he_api.cpp he_eos.h he_api.h
	$(CXX) $(CXXFLAGS) he_fit.cpp he_eos.cpp he_api.cpp -o $@

# таблица + полный набор проверок
test: all
	./he_table_build -o he_table2.bin --tol 1e-4 --selfcheck 2000
	./he_test he_table2.bin

clean:
	rm -f he_table_build he_calc he_test he_fit he_table2.bin
.PHONY: all coolprop test clean
