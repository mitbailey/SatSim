CXX = g++
CPPOBJS = src/main.o
COBJS = src/wmm_geomag.o
EDCXXFLAGS = -I ./ -I ./include/ -Wall -pthread $(CXXFLAGS)
EDLDFLAGS :=
TARGET = sim.out

all: $(COBJS) $(CPPOBJS)
	$(CXX) $(EDCXXFLAGS) $(COBJS) $(CPPOBJS) -o $(TARGET) $(EDLDFLAGS)
	sudo ./$(TARGET)

%.o: %.cpp
	$(CXX) $(EDCXXFLAGS) -o $@ -c $<

%.o: %.c
	$(CXX) $(EDCXXFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	$(RM) *.out
	$(RM) *.o
	$(RM) src/*.o