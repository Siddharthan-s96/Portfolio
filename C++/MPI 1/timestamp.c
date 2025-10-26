#include <stdio.h>
#include <sys/time.h>
#include <time.h>

int
main()
{
	char           buffer[26];
	struct tm*     tm_info;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	tm_info = localtime(&tv.tv_sec);

	strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
	printf("%s.%06ld\n", buffer, tv.tv_usec);

	return 0;
}
