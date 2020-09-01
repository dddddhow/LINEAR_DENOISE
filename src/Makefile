#makefile for c c++ fotran
divider="--------------------------------------------------------------------------------"
BinDir =../bin
LibDir =../lib
IncDir =../inc
ObjDir =../obj
DatDir =../file
BackDir=~/
Date=$(shell date +%F_time%H)
#--------------------------------------#
Projn   =LinearRadonTransform
#--------------------------------------#
CPP_CMPL_0=icpc
CPP_CMPL=$(CPP_CMPL_0) -DARMA_ALLOW_FAKE_GCC
C_CMPL  =gcc
F_CMPL  =gfortran
LINKER  =$(CPP_CMPL_0)

C_CMPL_FLAG   = $(addprefix -I, ${IncDir})
CPP_CMPL_FLAG = $(addprefix -I, ${IncDir})
F_CMPL_FLAG   =
LINKER_FLAG   = $(addprefix -L, ${LibDir} )  -lgfortran -larmadillo
LINKLIBS	  = -lm

CSrc   =$(wildcard *.c)
CppSrc =$(wildcard *.cpp)
F90Src =$(wildcard *.f90)
F77Src =$(wildcard *.f)

CObj0  =$(patsubst %.c,%.o,${CSrc})
CppObj0=$(patsubst %.cpp,%.o,${CppSrc})
F90Obj0=$(patsubst %.f90,%.o,${F90Src})
F77Obj0=$(patsubst %.f,%.o,${F77Src})

CObj   =$(addprefix $(ObjDir)/,$(CObj0))
CppObj =$(addprefix $(ObjDir)/,$(CppObj0))
F90Obj =$(addprefix $(ObjDir)/,$(F90Obj0))
F77Obj =$(addprefix $(ObjDir)/,$(F77Obj0))

OBJ    =$(CObj) $(CppObj) $(F90Obj) $(F77Obj)
EXE    =$(BinDir)/$(Projn).e
BACK   =$(BackDir)/$(Projn)_back_$(Date).tgz
Back_Content=$(shell basename $(shell pwd))
ex_tar ="--exclude=*.sw? --exclude=*.log --exclude=*.[oe]"

all:
	mkdir -p $(ObjDir)
	make -j 8 $(EXE)

$(EXE):$(OBJ)
	@echo $(divider)
	@echo '			Linking ---> ' $(EXE)
	@echo $(divider)
	#$(LINKER) $^ -o $(EXE)  $(LINKLIBS) -fopenmp ${LINKER_FLAG}
	$(LINKER) ${OBJ} -o $(EXE)  $(LINKLIBS) -fopenmp ${LINKER_FLAG}
$(ObjDir):
	mkdir -p $(ObjDir)
$(CppObj):$(ObjDir)/%.o:%.cpp
	@echo $(divider)
	@echo '			Making '$<'  --> *.o'
	@echo $(divider)
	$(CPP_CMPL) -c  ${CPP_CMPL_FLAG} $< -o $@ $(CPP_CMPL_FLAG)
$(CObj):$(ObjDir)/%.o:%.c
	@echo $(divider)
	@echo '			Making '$<'  --> *.o'
	@echo $(divider)
	$(C_CMPL)   -c  ${C_CMPL_FLAG} $< -o $@ -fopenmp
$(F90Obj):$(ObjDir)/%.o:%.f90
	@echo $(divider)
	@echo '			Making *.f90  --> *.o'
	@echo $(divider)
	$(F_CMPL)  -c  $< -o $@   -fopenmp -cpp
$(F77Obj):$(ObjDir)/%.o:%.f
	@echo $(divider)
	@echo '			Making *.f  --> *.o'
	@echo $(divider) -fopenmp -cpp
	$(F_CMPL)  -c  $< -o $@

backup:
	@#make clean;
	mkdir -p $(BackDir)
	cd .. && tar $(ex_tar) -czvf $(BACK) $(Back_Content)
	ls -hl $(BACK)
b:backup

clean:
	@#CREAT EMPTY files
	@touch $(EXE)
	@touch $(OBJ)
	rm  $(EXE) $(OBJ)
	@pwd; echo $(Date)
