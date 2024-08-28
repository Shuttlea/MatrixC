CC = gcc -Wall -Werror -Wextra -std=c11
CHECK_FLAGS = $(shell pkg-config --cflags --libs check)
COVERAGE = -fprofile-arcs -ftest-coverage

all: format gcov_report

test: s21_matrix.a s21_matrix_tests.c s21_matrix.c s21_matrix.h
	$(CC) $(COVERAGE) $(CHECK_FLAGS)  -c *.c -g
	$(CC) $(COVERAGE) -o s21_matrix_tests s21_matrix_tests.o s21_matrix.o s21_matrix.a $(CHECK_FLAGS) -lm
	./s21_matrix_tests

gcov_report: test
	gcovr --gcov-executable gcov-13 -r . --html --html-details -o coverage_report.html
	

s21_matrix.a: s21_matrix.c s21_matrix.h
	$(CC) -c s21_matrix.c
	ar crs s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a
	rm -rf *.o

format:
	clang-format -i *.h *.c

clean:
	rm -rf *.o *.a *.gcno *.gch *.gcda *.gcov *.css *.html *.out s21_matrix_tests