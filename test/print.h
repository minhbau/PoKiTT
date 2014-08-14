/*
 * print.h
 *
 *  Created on: May 20, 2014
 *      Author: nate
 */

#ifndef PRINT_H_
#define PRINT_H_

#include<vector>
#include<string>
#include<iostream>
#include<cmath>
#include<numeric>

void find(double place){
#ifdef FIND
  std::cout<<place<<std::endl;
#endif
}

void print(std::string label, std::vector<int> vector){
  int index=0;
  for( std::vector<int>::iterator it = vector.begin(); it!=vector.end(); ++it, ++index){
    std::cout<<label<<" "<<index<<" = "<<*it<<std::endl;
  }
  std::cout<<std::endl;
}

void print(std::string label, double value, int index=-1){
  if(index==-1){
    std::cout<<label<<" = "<<value<<std::endl;
  }
  else{
    std::cout<<label<<" "<<index<<" = "<<value<<std::endl;
  }
}

void print(std::string label, std::vector<double> vector){
  int index=0;
  for( std::vector<double>::iterator it = vector.begin(); it!=vector.end(); ++it, ++index){
    std::cout<<label<<" "<<index<<" = "<<*it<<std::endl;
  }
  std::cout<<std::endl;
}



#endif /* PRINT_H_ */
