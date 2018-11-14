#makefileのテンプレート
#汎用的に使えそうなので書いておく

CC = gcc
CFLAGS = -g -MMD -MP -Wall -Wextra -Winit-self -Who-missing-field-initializers
ifeq "$(shell getconf LONG_BIT)" "64"
	LDFLAGS = 
else
	LDFLAGS = 
endif
LIBS = -lm -lpthread 
INCLUDE = -I./include
TARGET = ./bin/$(shell basename `readlink -f .`)
SRCDIR = ./source
ifeq "$(strip $(SRCDIR))" ""
	SRCDIR = .
endif
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJDIR = ./obj
ifeq "$(strip $(SRCDIR))" ""
	OBJDIR = .
endif
OBJECTS = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.c=.o)))
DEPENDS = $(OBJECTS:.o=.d)

$(TARGET): $(OBJECTS) $(LIBS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	-mkdir -p $(OBJDIR)
	$(CC) $(CLAGS) $(INCLUDE) -o $@ -c $<

all: clean $(TARGET)

clean:
	-rm -f $(OBJECTS) $(DEPENDS) $(TARGET)

-include $(DEPENDS)
