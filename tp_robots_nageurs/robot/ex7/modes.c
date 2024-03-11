#include "config.h"
#include "modes.h"
#include "robot.h"
#include "module.h"
#include "registers.h"
#include "hardware.h"
#include "can.h"

volatile float freq  =  1.0;   // Hz
volatile float ampl_set  = 0.0;
volatile float ampl;
volatile float phase =  1.5;

const float AMPLITUDE_DECREASE_RATE = 0.1;

#define NUM_POS 5
const uint8_t addresses[NUM_POS] = {5,26,24,22,25}; // ToDo check addresses, to insert starting from tail


/* Register callback function, handles some new registers on the radio.
 * All these registers are of course completely useless, but it demonstrates how
 * to implement a register callback function, and what it can do.
 */
static int8_t register_handler(uint8_t operation, uint8_t address, RadioData* radio_data)
{
  switch (operation)
  {
    case ROP_WRITE_8:
      if (address == 10) {
        freq = DECODE_PARAM_8(radio_data->byte,(0),(2));

        if(freq > 2.0f) freq = 2.0f;
        
        return TRUE;
      }
      else if (address == 11) {
        ampl_set = DECODE_PARAM_8(radio_data->byte,(0),(60));

        if(ampl_set < 0.0f) ampl_set = 0.0f;
        else if(ampl_set > 60.0f) ampl_set = 60.0f;

        return TRUE;
      } else if (address == 12) {
        phase = DECODE_PARAM_8(radio_data->byte,(0.5),(1.5));

        if(phase < 0.5f) phase = 0.5f;
        else if(phase > 1.5f) phase = 1.5f;

        return TRUE;
      }
  }
  return FALSE;
}



void traveling_wave_demo_mode()
{
  uint32_t dt, cycletimer;
  float my_time, delta_t, l;
  int8_t l_rounded;

  cycletimer = getSysTICs();
  my_time = 0;

  for(uint8_t i = 0; i < NUM_POS; i++){
    init_body_module(addresses[i]);
    set_reg_value_dw(addresses[i], MREG32_LED, 0);
    start_pid(addresses[i]);
  }  

  do {
    // Calculates the delta_t in seconds and adds it to the current time
    dt = getElapsedSysTICs(cycletimer);
    cycletimer = getSysTICs();
    delta_t = (float) dt / sysTICSperSEC;


    my_time += delta_t;

    for (uint8_t i = 0; i < NUM_POS; i++) {
      ampl = (1 - i * AMPLITUDE_DECREASE_RATE) > 0 ? ampl_set * (1 - i * AMPLITUDE_DECREASE_RATE > 0) : 0;
      
      l = ampl * sin(M_TWOPI * (freq * my_time + i * phase / NUM_POS));
      l_rounded = (int8_t) l;

      bus_set(addresses[i], MREG_SETPOINT, DEG_TO_OUTPUT_BODY(l_rounded));
    }

    
    // Make sure there is some delay, so that the timer output is not zero
    pause(ONE_MS);

  } while (reg8_table[REG8_MODE] == IMODE_TRAVELING_WAVE_DEMO);

  for (uint8_t i = 0; i < NUM_POS; i++)
    bus_set(addresses[i], MREG_SETPOINT, DEG_TO_OUTPUT_BODY(0.0));
  
  pause(ONE_SEC);
  
  // for (uint8_t i = 0; i < NUM_POS; i++)
  //   bus_set(addresses[i], MREG_MODE, MODE_IDLE);
  
  // Back to the "normal" green
  set_color(2);
}

void main_mode_loop()
{
  reg8_table[REG8_MODE] = IMODE_IDLE;
  
  // Registers the register handler callback function
  radio_add_reg_callback(register_handler);

  while (1)
  {
    switch(reg8_table[REG8_MODE])
    {
      case IMODE_IDLE:
        break;
      case IMODE_TRAVELING_WAVE_DEMO:
        traveling_wave_demo_mode();
        break;
      default:
        reg8_table[REG8_MODE] = IMODE_IDLE;
    }
  }
}
