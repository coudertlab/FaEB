all: src/main.cpp
	c++ -std=c++11 -I./include -O3 src/barrier.cpp -o faeb
clean:
	rm -f faeb
