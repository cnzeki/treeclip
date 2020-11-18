#include "ticker.h"

static Ticker s_ticker;

Ticker& GetTicker() {
	return s_ticker;
}
