#pragma once

#include <QString>
#include <QMap>
#include <QElapsedTimer>
#include <stdio.h>
#include <time.h>

class Stats
{

public:
	QString label;

	QElapsedTimer timer;
	int time;
	bool noTime;
	double value;
	bool ended;

public:
	Stats();
	Stats(QString stat_label);
	Stats(QString stat_label, double newValue);

	Stats(const Stats& from);
	Stats& operator= (const Stats& from);

	void start();
	void end();

	void addValue(double newValue);

	void print();

	double getValue();
	QString getLabel();
	bool isTimeBased();
};

class StatCollection
{
public:
	StatCollection(QString title);

	QMap<QString, Stats> myStats;
	Stats total;

	void start(QString stat_lable);
	void end(QString stat_lable);

	void print();
};
