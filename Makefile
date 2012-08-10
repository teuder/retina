# コンパイラとオプション
# #includeのパスを追加。頭に"-I"をつける。例 := -I../header
INCLUDES :=

CXX :=g++-4.6



CXXFLAGS_D := -g -Wall -std=c++0x -DDEBUG  -frounding-math $(INCLUDES)
CXXFLAGS_R := -Wall -std=c++0x -DNDEBUG -frounding-math -O3 $(INCLUDES) 


#CXX :=clang++
#CXXFLAGS_D := -Wall -std=c++0x -DDEBUG  -frounding-math -stdlib=libc++ $(INCLUDES)
#CXXFLAGS_R := -Wall -std=c++0x -DNDEBUG -DCGAL_CFG_NO_TR1_ARRAY -DCGAL_CFG_NO_TR1_TUPLE -frounding-math -stdlib=libc++ -O3 $(INCLUDES) 



# ソースファイルの指定
SRCS := $(wildcard *.cpp)
OBJS_D := $(SRCS:.cpp=_D.o) #for release
OBJS_R := $(SRCS:.cpp=_R.o) #for debug


# 最終産物の名前
PROGRAM := retina

# 生成物を置くディレクトリの指定
OBJD := ./
PARD := ~/bin/
OBJS_D := $(addprefix $(OBJD),$(OBJS_D))
OBJS_R := $(addprefix $(OBJD),$(OBJS_R))

PROGRAM := $(PARD)$(PROGRAM)
OBJDESC := $(subst /,\/,$(OBJD))




#リンクするライブラリ
LIBS := -lgflags -lCGAL -lgmp -lboost_thread -lglut -lGLU -lGL

# コンパイルせずリンクして使うだけの既製オブジェクトを指定
EXOBJS := ./SFMT/SFMT.o



# デフォルトターゲット（ただmakeしたとき実行される)
.DEFAULT_GOAL := release

# 依存関係
DEPS := Dependfile		# ファイル名
include $(DEPS)		# 読み込み

# 依存関係自動生成用ターゲット（SRCSを読んでDEPSに書き出す）
.PHONY: depend
depend: 
	$(CXX) -MM $(INCLUDES) $(SRCS) | sed 's/^\(.*:\)/$(OBJDESC)\1/' > $(DEPS);
	$(CXX) -MM $(INCLUDES) $(SRCS) | sed 's/^\(.*\)\(.o:\)/$(OBJDESC)\1_R\2/' >> $(DEPS);
	$(CXX) -MM $(INCLUDES) $(SRCS) | sed 's/^\(.*\)\(.o:\)/$(OBJDESC)\1_D\2/' >> $(DEPS);



# 最終産物コンパイル用ターゲット(make release)
.PHONY: debug
debug: $(OBJS_D) $(EXOBJS)
	$(CXX) $(CXXFLAGS_D) -o $(PROGRAM) $^ $(LIBS)

.PHONY: release
release: $(OBJS_R) $(EXOBJS)
	$(CXX) $(CXXFLAGS_R) -o $(PROGRAM) $^ $(LIBS)



# 生成物削除用ターゲット（make clean）
.PHONY: clean
clean:
	rm -rf $(OBJS_R) $(OBJS_D) $(PROGRAM)

# パターンルール（*.cppから*.oを作れ）
$(OBJD)%_D.o: %.cpp
	$(CXX) $(CXXFLAGS_D) -o $@ -c $<

$(OBJD)%_R.o: %.cpp
	$(CXX) $(CXXFLAGS_R) -o $@ -c $<

.PHONY: SFMT
SFMT:
	g++ -c -O3 -msse2 -DHAVE_SSE2=1 -DMEXP=19937 ./SFMT/SFMT.c -o ./SFMT/SFMT.o





