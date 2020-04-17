import pytest
from merge_metadata import *

def test_format_dates():
	dates = ["2020-01-1", "2020-01-01", \
		"2020-01", "2020", "Jan 2020", \
		"Jan 01 2020", "Jan-2020", \
		"Jan-01-2020"]
	formatted = format_dates(dates)
	assert formatted == ["2020-01-01", "2020-01-01", \
		"2020-01", "2020", "2020-01", \
		"2020-01-01", "2020-01", \
		"2020-01-01"]

def test_parse_location():
	locations = ["Europe / United Kingdom / England", \
		"United Kingdom / England", \
		"Africa / Algeria / Boufarik", \
		"Europe / Italy / Abruzzo / Teramo", \
		"Europe / France / Hauts de France / Compi√®gne", \
		"Asia / India", \
		"Oceania / Australia / Victoria", \
		"Africa / Democratic Republic of the Congo / Kinshasa", \
		"Asia / Hong Kong", \
		"", \
		"USA:NC", \
		"China:Zhejiang, Hangzhou", \
		"USA:CA, San Diego County"]
	country, region = parse_location(locations)
	assert country == ["United Kingdom", \
		"United Kingdom", \
		"Algeria", \
		"Italy", \
		"France", \
		"India", \
		"Australia", \
		"Democratic Republic of the Congo", \
		"China", \
		"", \
		"USA", \
		"China", \
		"USA"]
	assert region == ["England", \
		"England", \
		"Boufarik", \
		"Abruzzo", \
		"Hauts de France", \
		"", \
		"Victoria", \
		"Kinshasa", \
		"Hong Kong", \
		"", \
		"NC", \
		"Zhejiang, Hangzhou", \
		"CA, San Diego County"]