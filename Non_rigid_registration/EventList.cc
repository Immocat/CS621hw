#include "EventList.hh"
void EventList::init(size_t frame_count){
  events_at_frame.clear();
  events_at_frame.assign(frame_count, std::vector<Event>());
}