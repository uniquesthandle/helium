CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++20 -Wall -Wextra
CORE      = he_api.cpp he_eos.cpp he_table.cpp

all: he_table_build he_calc he_test he_fit he_check_abinitio

he_table_build: he_table_build.cpp he_ref.cpp $(CORE) he_api.h he_eos.h he_table.h he_ref.h he_console.h
	$(CXX) $(CXXFLAGS) he_table_build.cpp he_ref.cpp $(CORE) -o $@

# Построитель с эталоном CoolProp (пути укажите свои):
#   make coolprop CP_INC="/e/misc/my_coolprop_project/externals/CoolProp/include" CP_LIB="/e/misc/my_coolprop_project/build/externals/CoolProp"

CP_INC = /c/code/my_coolprop_project/externals/CoolProp/include -I/c/code/my_coolprop_project/build/_deps/fmt-src/include
CP_LIB = /c/code/my_coolprop_project/build/CoolProp/Release

coolprop: he_table_build.cpp he_ref.cpp $(CORE) he_console.h
	$(CXX) $(CXXFLAGS) -DUSE_COOLPROP -I$(CP_INC) he_table_build.cpp he_ref.cpp $(CORE) -L$(CP_LIB) -lCoolProp -o he_table_build

he_calc: he_calc.cpp $(CORE) he_api.h he_eos.h he_table.h
	$(CXX) $(CXXFLAGS) he_calc.cpp $(CORE) -o $@

he_test: he_test.cpp $(CORE) he_api.h he_eos.h he_table.h
	$(CXX) $(CXXFLAGS) he_test.cpp $(CORE) -o $@

he_fit: he_fit.cpp he_eos.cpp he_api.cpp he_eos.h he_api.h
	$(CXX) $(CXXFLAGS) he_fit.cpp he_eos.cpp he_api.cpp -o $@

he_check_abinitio: he_check_abinitio.cpp $(CORE) he_api.h he_eos.h he_table.h
	$(CXX) $(CXXFLAGS) he_check_abinitio.cpp $(CORE) -o $@

he_gas_test: he_gas_test.cpp he_gas.cpp he_api.cpp he_table.cpp he_eos.cpp he_console.h
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^

# проверка по ab initio-значениям (Bich, Hellmann, Vogel 2007): встроенные коэффициенты (ожидаемо не проходят),
# затем привязка к ab initio и повторная проверка
check-abinitio: he_check_abinitio he_fit
	-./he_check_abinitio data/he4_abinitio_bich2007.csv
	./he_fit - --abinitio data/he4_abinitio_bich2007.csv -o he_coef_abinitio.txt
	./he_check_abinitio data/he4_abinitio_bich2007.csv --model analytic:he_coef_abinitio.txt

# таблица + полный набор проверок
test: all
	./he_table_build -o he_table2.bin --tol 1e-4 --selfcheck 2000
	./he_test he_table2.bin

clean:
	rm -f he_table_build he_calc he_test he_fit he_check_abinitio he_table2.bin he_coef_abinitio.txt
.PHONY: all coolprop test check-abinitio clean
