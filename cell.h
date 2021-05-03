/* 
 * File:   cell.h
 * Author: Peter Kováč
 *
 * Created on April 26, 2016, 9:17 PM
 */

#ifndef _CELL_H_
#define _CELL_H_

class Cell
{
public:
	///Constructor of Cell
	Cell(): volume(0) { }
  
	/// Return cell's current volume
	int Volume() const {
		return volume;
	}
	
	/// Decrement cell's volume by one unit
	int IncrementVolume() {
		return ++volume;
	}
  
	/// Increment cell's volume by one unit
	int DecrementVolume() {
		return --volume;
	}
  
private:
	/// current volume of cell
	int volume;
	
};

#endif	/* _CELL_H_ */

