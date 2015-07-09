basedirs      = src exe
macrodirs     = macros
testdirs      = test

.PHONY: base macros clean 

all: base

base: 
	@for dir in $(basedirs); do (cd $$dir; make); done

macros: 
	@for dir in $(macrodirs); do (cd $$dir; make); done

test: 
	@for dir in $(testdirs); do (cd $$dir; make); done

clean:
	@rm -f *~
	@for dir in $(basedirs); do (cd $$dir; make clean; rm -f *~); done
	@for dir in $(macrodirs); do (cd $$dir; make clean; rm -f *~); done
	@for dir in $(testdirs); do (cd $$dir; make clean; rm -f *~); done

