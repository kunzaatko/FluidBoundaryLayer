default: clean build

LIBS:="-lm"

CFLAGS := "-std=c99 " + "-g " + "-Wall " + "-Wextra " + \
"-pedantic " +"-Werror " +"-Wmissing-declarations " + "-DUNITY_SUPPORT_64 "

ASANFLAGS  := "-fsanitize=address " + "-fno-common " + "-fno-omit-frame-pointer "
alias b:=build

build: build-program build-diagrams build-report build-plots

build-program: 
    @echo Building program...
    cc {{ CFLAGS }} src/*.c -o program {{ LIBS }}

build-diagrams:
    @echo Building diagrams...
    cd ./report/diagrams/ && tectonic -X build # LaTeX build

    # Converting to PNG
    convert -density 720 -background transparent -flatten ./report/diagrams/build/wedge/wedge.pdf  static/wedge.png
    convert -density 720 -background transparent -flatten ./report/diagrams/build/flow_profile/flow_profile.pdf
    static/flow_profile.png

build-report:
    @echo Building report...
    cd ./report/report/ && tectonic -X build # LaTeX build

build-plots:
    @echo Building report...
    cd ./report/plots/ && tectonic -X build # LaTeX build

alias t:=test
test: testprog
    @echo Running tests...
    ./tests.out

testprog:
    @echo Compiling...
    cc {{ CFLAGS }} src/util.c src/blais.c test/vendor/unity.c test/*.c -o tests.out {{LIBS}}

alias c:=clean

clean: clean-program clean-diagrams clean-report clean-plots
clean-program:
    @echo Cleaning program...
    rm -rf *.o *.out *.out.dSYM program

clean-diagrams:
    @echo Cleaning diagrams...
    rm -f static/*.png
    rm -rf ./report/diagrams/build/*/*.pdf

clean-report:
    @echo Cleaning report...
    rm -rf ./report/protocol/build/*/*.pdf

clean-plots:
    @echo Cleaning plots...
    rm -rf ./report/plots/build/*/*.pdf
