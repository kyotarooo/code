#コンパイラ指定
CC = gcc
#コンパイルオプションの指定
CFLAGS = -O3 
#実行ファイル名
TARGETS = grid_vtk merge_grid section_vtk

all: $(TARGETS)
#実行ファイルに名前をつける時 -o
grid_vtk: grid_vtk.c
	$(CC) $(CFLAGS) -o grid_vtk grid_vtk.c
#:の左　生成物の名前
#:の右　依存先
merge_grid: merge_grid.c
	$(CC) $(CFLAGS) -o merge_grid merge_grid.c

section_vtk: section_vtk.c
	$(CC) $(CFLAGS) -o section_vtk section_vtk.c -lm

clean:
	rm -f $(TARGETS) *.o
#clean-allはcleanに連動して動く
clean-all: clean
	rm -f *.vtk 




