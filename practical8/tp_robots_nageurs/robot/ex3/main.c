#include "hardware.h"
#include "module.h"
#include "robot.h"
#include "registers.h"

#define NUM_POS 4
const uint8_t addresses[NUM_POS] = {72,73,74,21};
int8_t positions[NUM_POS];

static int8_t register_handler(uint8_t operation, uint8_t address, RadioData* radio_data)
{
  if(operation == ROP_READ_MB && address == 2)
  {
    radio_data->multibyte.size = NUM_POS;
    for (uint8_t i = 0; i < NUM_POS; i++) 
      radio_data->multibyte.data[i] = positions[i];
    return TRUE;
  }
  return FALSE;
}

int main(void)
{
  hardware_init();
  
  // Registers the register handler callback function
  radio_add_reg_callback(register_handler);
  
  // Changes the color of the led (red) to show the boot
  set_color_i(4, 0);

  // Initialises the body module with the specified address   
  init_body_module(addresses[0]);
  init_body_module(addresses[3]);
  init_limb_module(addresses[1]);
  init_limb_module(addresses[2]);
  
  while (1) 
  {
    // Get the position of each joint
    for(uint8_t i = 0; i < NUM_POS; i++)
      positions[i] = bus_get(addresses[i], MREG_POSITION);
  }

  return 0;
}
