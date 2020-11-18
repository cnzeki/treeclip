#ifndef TICKER_H
#define TICKER_H

#include <chrono> 
#include <string>
#include <vector>
#include <map>
using namespace std::chrono;

#define TTAG_FIRST  "__FIRST__"
#define TTAG_LAST   "__LAST__"
#define TTAG_NOW    "__NOW__"

#define system_clock high_resolution_clock

class Ticker {
    // ticks
	std::map< std::string, system_clock::time_point > ticks;
	std::vector< std::string > tickTags;
	// tick subjects
	std::vector< std::pair<std::string, float> > records;

	void gettime(system_clock::time_point &t)
	{
		t = system_clock::now();
	}

public:

	Ticker() { }

	void reset() {
		records.clear();
		ticks.clear();
	}

	float getTime(std::string tickTag)
	{
		for (int i = 0; i < records.size(); i++)
		{
			if (tickTag == records[i].first)
				return records[i].second;
		}

		return 0.;
	}

	void tick(const char* tickTag) {
		system_clock::time_point now;
		gettime(now);
		ticks.insert(std::pair<std::string, system_clock::time_point>(
			std::string(tickTag), now));
		
		tickTags.push_back(std::string(tickTag));
	}

	void record(const char* recTag, const char* tickStart = TTAG_LAST,
		const char* tickStop=TTAG_NOW) {
		std::string start(tickStart);
		if (start == TTAG_LAST) {
			start = tickTags[tickTags.size()-1];
		}
		else if (start == TTAG_FIRST)
		{
			start = tickTags[0];
		}

		std::string stop(tickStop);
		if (stop == TTAG_NOW) {
			tick(recTag);
			stop = std::string(recTag);
		}
		
		std::map< std::string, system_clock::time_point >::iterator iterStart,
			iterStop;
		iterStart = ticks.find(start);
		if (iterStart == ticks.end())
		{
			printf("Start tick not found: %s\n", tickStart);
		}
		iterStop = ticks.find(stop);
		if (iterStop == ticks.end())
		{
			printf("Stop tick not found: %s\n", tickStop);
		}
		// record time
		system_clock::time_point timeStart = iterStart->second;
		system_clock::time_point timeStop = iterStop->second;
		system_clock::duration t = timeStop - timeStart;
		std::pair<std::string, float> record(std::string(recTag), toMSecs(t));
		records.push_back(record);
	}

	void report(const char* title = NULL) {
		if (title) {
			printf("Timing report for %s:\n", title);
		}
		else
		{
			printf("Records:\n");
		}

		for (int i = 0; i < records.size(); i++)
		{
			printf("  %20s: %10.3f\n", records[i].first.c_str(), records[i].second);
		}
		printf("\n");
	}

	void reportTicks() {
		printf("Ticks  :\n");
		std::map< std::string, system_clock::time_point >::iterator iterStart, iter;
		iterStart = ticks.find(tickTags[0]);
		system_clock::time_point timeStart = iterStart->second;
		for (int i = 0; i < tickTags.size(); i++)
		{
			iter = ticks.find(tickTags[i]);
			system_clock::time_point timeStop = iter->second;
			system_clock::duration t = timeStop - timeStart;
			printf("  %20s: %10.3f\n", iter->first.c_str(), toMSecs(t));
		}
	}

	float toSecs(system_clock::duration& d) { return toUSecs(d) * 1e-6; }
	float toMSecs(system_clock::duration& d) { return toUSecs(d) * 1e-3;}
	float toUSecs(system_clock::duration& d) { return float(duration_cast<microseconds>(d).count()); }
};

Ticker& GetTicker();

#endif
