CC      := clang
CCFLAGS := 
LDFLAGS :=

TARGETS:= fft-rs
OBJ    := RSErasureCode.o
# DEPS   := RSErasureCode.h

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJ)

$(OBJ): %.o : %.c 
	$(CC) -c -o $@ $< $(CCFLAGS)

$(TARGETS): % : $(OBJ)
	$(CC) -o $@ $(LIBS) $^ $(CCFLAGS) $(LDFLAGS)

run: $(TARGETS)
	valgrind ./fft-rs
