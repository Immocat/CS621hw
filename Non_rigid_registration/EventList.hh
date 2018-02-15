// Chapter 4.4 of "Tracking Surfaces with Evolving Topology"
#pragma once
#include <vector>
class Event {
 public:
  enum EType{
    SPLIT_EDGE,
    SPLIT_TRIANGLE,
    COLLAPSE_EDGE,
    NUM_OF_EVENT_TYPE
  };
  Event() {}
  public:
  //public data
  EType etype;
  unsigned int start[3];
  unsigned int end;
  double alpha[3];

};
class EventList {
 public:
  EventList() {}
  void init(size_t frames);
  std::vector<std::vector<Event>> events_at_frame;
};