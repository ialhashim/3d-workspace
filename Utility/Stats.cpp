#include "Stats.h"

Stats::Stats( QString stat_label )
{
	this->label = stat_label;
	time = 0;
	noTime = false;
	ended = true;

	// Auto start
	start();
}

Stats::Stats()
{
	this->label = "NOTHING";
	time = 0;
	noTime = false;
	ended = true;
}

Stats::Stats( QString stat_label, double newValue )
{
	this->label = stat_label;
	time = 0;
	noTime = true;
	value = newValue;
	ended = true;
}

Stats::Stats( const Stats& from )
{
	this->label = from.label;
	this->time = from.time;
	this->value = from.value;
	this->noTime = from.noTime;
	this->ended = from.ended;
}

Stats& Stats::operator=( const Stats& from )
{
	this->label = from.label;
	this->time = from.time;
	this->value = from.value;
	this->noTime = from.noTime;
	this->ended = from.ended;

	return *this;
}

void Stats::start()
{
	this->timer.start();

	this->ended = false;
}

void Stats::end()
{
	if(!ended)
	{
		this->time = (int)timer.elapsed();
		ended = true;
	}
}

void Stats::print()
{
	end();

	if(noTime)
		printf(" :: %s \t =   %f  \n", label.toStdString().c_str(), this->value);
	else
		printf(" :: %s \t time =  %d   ms\n", label.toStdString().c_str(), (int)time);
}

void Stats::addValue( double newValue )
{
	this->value = newValue;
	this->noTime = true;
}

double Stats::getValue()
{
	if(noTime)
		return value;
	else
		return (double) time;
}

QString Stats::getLabel()
{
	return label;
}

bool Stats::isTimeBased()
{
	return !noTime;
}

// For collection of statistics
StatCollection::StatCollection( QString title )
{
	total = Stats(title);
	total.start();
}
void StatCollection::start( QString stat_label )
{
	myStats[stat_label] = Stats(stat_label);
	myStats[stat_label].start();
}

void StatCollection::end( QString stat_label )
{
	myStats[stat_label].end();
}

void StatCollection::print()
{
	total.end();

	printf(":: %s ", total.label.toStdString().c_str());

	foreach (Stats s, myStats)
	{
		printf("+ %s (%d ms)", s.label.toStdString().c_str(), (int)s.time);
	}

	printf(" = %d ms\n", (int)total.time);
}
