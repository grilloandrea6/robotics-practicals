#include "hardware.h"
#include "module.h"
#include "robot.h"
#include "registers.h"
#include "can.h"

#define NUM_POS 4
const uint8_t addresses[NUM_POS] = {72,73,74,21};



int main(void)
{
  hardware_init();
  
  // Initialises the body module with the specified address (but do not start
  // the PD controller)
  
  // ToDo set addresses and check number of modules
  for(uint8_t i = 0; i < NUM_POS; i++) {
    init_body_module(addresses[i]);
    set_reg_value_dw(addresses[i], MREG32_LED, 0);
  }  

  return 0;
}
