#include "hardware.h"
#include "module.h"
#include "robot.h"
#include "registers.h"

#define NUM_POS 4
const uint8_t addresses[NUM_POS] = {72,73,74,21};
int8_t positions[NUM_POS];


/* Register callback function, handles some new registers on the radio.
 * All these registers are of course completely useless, but it demonstrates how
 * to implement a register callback function, and what it can do.
 */
static int8_t register_handler(uint8_t operation, uint8_t address, RadioData* radio_data)
{
  uint8_t i;
  
  switch (operation)
  {
    case ROP_READ_MB:
      if (address == 2) {
        radio_data->multibyte.size = NUM_POS;
        for (i = 0; i < NUM_POS; i++) {
          radio_data->multibyte.data[i] = positions[i];
        }
        return TRUE;
      }
      break;
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

  // Initialises the body module with the specified address (but do not start
  // the PD controller)
  
  init_body_module(addresses[0]);
  init_body_module(addresses[3]);
  init_limb_module(addresses[1]);
  init_limb_module(addresses[2]);
  
  
  // And then... do this
  while (1) {
    for(uint8_t i = 0; i < NUM_POS; i++) {
      positions[i] = bus_get(addresses[i], MREG_POSITION);
    }
    // pause(TEN_MS); // not to flood canbus
  }

  return 0;
}
