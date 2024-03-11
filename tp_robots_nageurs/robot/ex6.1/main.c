#include "hardware.h"
#include "module.h"
#include "robot.h"
#include "registers.h"
#include "can.h"

#define NUM_POS 5
const uint8_t addresses[NUM_POS] = {25,22,24,26,5};

int main(void)
{
  hardware_init();
  
  // Initialises all the body modules with the specified address (but do not start
  // the PD controllers)
  for(uint8_t i = 0; i < NUM_POS; i++) {
    init_body_module(addresses[i]);
    set_reg_value_dw(addresses[i], MREG32_LED, 0);
  }  

  return 0;
}
