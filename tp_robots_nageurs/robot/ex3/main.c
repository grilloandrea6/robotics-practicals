#include "hardware.h"
#include "module.h"
#include "robot.h"
#include "registers.h"

const uint8_t MOTOR_ADDR = 21;

int8_t positions[3];


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
        radio_data->multibyte.size = 3;
        for (i = 0; i < 3; i++) {
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
  
  // Changes the color of the led (red) to show the boot
  set_color_i(4, 0);

  // Initialises the body module with the specified address (but do not start
  // the PD controller)
  init_body_module(MOTOR_ADDR);
  init_body_module(MOTOR_ADDR + 1); // is it needed?
  init_body_module(MOTOR_ADDR + 2); // is it needed?
  
  // And then... do this
  while (1) {
    positions[0] = bus_get(MOTOR_ADDR, MREG_POSITION);
    positions[1] = bus_get(MOTOR_ADDR + 1, MREG_POSITION);
    positions[2] = bus_get(MOTOR_ADDR + 2, MREG_POSITION);
    pause(TEN_MS); // not to flood canbus
  }

  return 0;
}
