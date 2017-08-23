/* 
* icl::protoip
* Authors: asuardi <https://github.com/asuardi>, bulatkhusainov <https://github.com/bulatkhusainov>
* Date: November - 2014
*/


#include <stdio.h>
#include <string.h>

#include "lwip/err.h"
#include "lwip/tcp.h"
#include "lwip/udp.h"
#ifdef __arm__
#include "xil_printf.h"
#endif
#include <stdint.h>

#include "xil_cache.h"
#include "xparameters.h"
#include "xparameters_ps.h"	/* defines XPAR values */
#include "xstatus.h"
#include "xil_io.h"
#include "xscutimer.h"
#include "platform_config.h"

#include "xfoo.h"
#include "FPGAserver.h"
#include "soc_user.h"

#define TIMER_DEVICE_ID		XPAR_XSCUTIMER_0_DEVICE_ID
#define TIMER_LOAD_VALUE	0xFFFFFFE
#define TIMER_RES_DIVIDER 	40
#define NSECS_PER_SEC 		667000000/2
#define EE_TICKS_PER_SEC 	(NSECS_PER_SEC / TIMER_RES_DIVIDER)

typedef uint32_t           Xint32;     /**< signed 32-bit */


//define as global variables
XFoo xcore;
//XFoo_Config config;
XScuTimer Timer;

unsigned int CntValue1 = 0;
unsigned int CntValue2 = 0;

Xint32 *pl_vec_in_in_ptr_ddr = (Xint32 *)pl_vec_in_IN_DEFINED_MEM_ADDRESS;
Xint32 *pl_mat_in_in_ptr_ddr = (Xint32 *)pl_mat_in_IN_DEFINED_MEM_ADDRESS;
Xint32 *pl_vec_out_out_ptr_ddr = (Xint32 *)pl_vec_out_OUT_DEFINED_MEM_ADDRESS;

/* Variables used for handling the packet */
uint8_t *payload_ptr;		/* Payload pointer */
Xint32   payload_temp;		/* 32-bit interpretation of payload */
uint8_t *payload_temp_char = (uint8_t *)&payload_temp;  /* Char interpretation of payload */

Xint32 inputvec[ETH_PACKET_LENGTH];
Xint32 outvec[ETH_PACKET_LENGTH];
int32_t inputvec_fix[ETH_PACKET_LENGTH];
int32_t outvec_fix[ETH_PACKET_LENGTH];

float soc_state0_in[SOC_STATE0_IN_VECTOR_LENGTH];
float soc_primal0_in[SOC_PRIMAL0_IN_VECTOR_LENGTH];
float soc_dual0_in[SOC_DUAL0_IN_VECTOR_LENGTH];
float soc_tol_iterates_in[SOC_TOL_ITERATES_IN_VECTOR_LENGTH];

float soc_primal_out[SOC_PRIMAL_OUT_VECTOR_LENGTH];
float soc_dual_out[SOC_DUAL_OUT_VECTOR_LENGTH];
float soc_aux_primal_out[SOC_AUX_PRIMAL_OUT_VECTOR_LENGTH];
float soc_iterates_out[SOC_ITERATES_OUT_VECTOR_LENGTH];

unsigned int write_offset;
int packet_internal_ID_previous;


int transfer_data() {
	return 0;
}


void print_app_header()
{
	xil_printf("\n\r");
	xil_printf("\n\r");
	xil_printf("---------------------------- ICL::PROTOIP-----------------------------------\n\r");
	xil_printf("------------------- (FPGA UPD/IP and TCP/IP server) ------------------------\n\r");
	xil_printf("----------------- asuardi <https://github.com/asuardi> ---------------------\n\r");
	xil_printf("-------------------------------- v1.0 --------------------------------------\n\r");
	xil_printf("\n\r");
	xil_printf("\n\r");
	if (TYPE_ETH==1)
		xil_printf("Starting TCP/IP server ...\n\r");
	if (TYPE_ETH==0)
		xil_printf("Starting UDP/IP server ...\n\r");
	xil_printf("\n\r");
	xil_printf("\n\r");
}


inline static void get_payload(){
	payload_temp_char[0] = *payload_ptr++;
	payload_temp_char[1] = *payload_ptr++;
	payload_temp_char[2] = *payload_ptr++;
	payload_temp_char[3] = *payload_ptr++;
}


unsigned int IpTimer(unsigned int TimerIntrId,unsigned short Mode)

{

    int                 Status;
    XScuTimer_Config    *ConfigPtr;
    volatile unsigned int     CntValue  = 0;
    XScuTimer           *TimerInstancePtr = &Timer;


    if (Mode == 0) {

      // Initialize the Private Timer so that it is ready to use

      ConfigPtr = XScuTimer_LookupConfig(TimerIntrId);

      Status = XScuTimer_CfgInitialize(TimerInstancePtr, ConfigPtr,
                     ConfigPtr->BaseAddr);

      if (Status != XST_SUCCESS) {
          return XST_FAILURE; }

      // Load the timer prescaler register.

      XScuTimer_SetPrescaler(TimerInstancePtr, TIMER_RES_DIVIDER);

      // Load the timer counter register.

      XScuTimer_LoadTimer(TimerInstancePtr, TIMER_LOAD_VALUE);

      // Start the timer counter and read start value

      XScuTimer_Start(TimerInstancePtr);
      CntValue = XScuTimer_GetCounterValue(TimerInstancePtr);

    }

    else {

       //  Read stop value and stop the timer counter

       CntValue = XScuTimer_GetCounterValue(TimerInstancePtr);
       XScuTimer_Stop(TimerInstancePtr);


    }

    return CntValue;

}


void udp_server_function(void *arg, struct udp_pcb *pcb,
		struct pbuf *p, struct ip_addr *addr, u16_t port){


	struct pbuf pnew;
	int k1;
	int k;
	int i;

	int packet_type;
	int packet_internal_ID;
	int packet_internal_ID_offset;


	float tmp_float;



	xcore.Bus_a_BaseAddress = 0x43c00000;
	xcore.IsReady = XIL_COMPONENT_IS_READY;



	// Only respond when the packet is the correct length
		if (p->len == (ETH_PACKET_LENGTH)*sizeof(int32_t)){


			/* Pick up a pointer to the payload */
			payload_ptr = (unsigned char *)p->payload;

			//Free the packet buffer
			pbuf_free(p);

			// Get the payload out
			for(k1=0;k1<ETH_PACKET_LENGTH;k1++){
				get_payload();
				inputvec[k1] = payload_temp;
			}

				//extract informations form the input packet
				tmp_float=*(float*)&inputvec[ETH_PACKET_LENGTH-2];
				packet_type=(int)tmp_float & 0x0000FFFF;
				packet_internal_ID=((int)tmp_float & 0xFFFF0000) >> 16; //if write packet_type, packet_num is the data vector ID

				tmp_float=*(float*)&inputvec[ETH_PACKET_LENGTH-1];
				packet_internal_ID_offset=(int)tmp_float; //if write packet_type, packet_num is the data vector ID

			if (DEBUG){
				printf("\n");
				printf("Received packet:\n");
				printf("packet_type=%x\n",packet_type);
				printf("packet_internal_ID=%x\n",packet_internal_ID);
				printf("packet_internal_ID_offset=%x\n",packet_internal_ID_offset);
			}


			if (packet_type==1) //reset IP
						{

							if (DEBUG)
								printf("Reset IP ...\n");

							if (XFoo_IsIdle(&xcore)==1)
							{
								if (DEBUG){
									printf("The core is ready to be used\n");
								}
							}
							else
								printf("ERROR: reprogram the FPGA !\n"); //should be added the IP reset procedure

							//Initialise output vector
							for(i=0;i<ETH_PACKET_LENGTH;i++)
							{
								tmp_float=1.0;
								outvec[i]=*(Xint32*)&tmp_float;
							}
							
							// Send back a dummy UDP packet
							pnew.next = NULL;
							pnew.payload = (unsigned char *)outvec;
							pnew.len = ETH_PACKET_LENGTH*sizeof(Xint32);
							pnew.type = PBUF_RAM;
							pnew.tot_len = pnew.len;
							pnew.ref = 1;
							pnew.flags = 0;

							udp_sendto(pcb, &pnew, addr, port);

						}
					else if (packet_type==2) //start SOC function
							{

								if (DEBUG)
									printf("Start SOC function ...\n");


								//Start timer
								CntValue1 = IpTimer(TIMER_DEVICE_ID,0);

								//Start SOC function
								soc_user(soc_state0_in,soc_primal0_in,soc_dual0_in,soc_tol_iterates_in,soc_primal_out,soc_dual_out,soc_aux_primal_out,soc_iterates_out);

								//Stop timer
								CntValue2 = IpTimer(TIMER_DEVICE_ID,1);

								if (DEBUG)
								{
									printf ("SOC Timer: beginning of the counter is : %d \n", CntValue1);
									printf ("SOC Timer: end of the counter is : %d \n", CntValue2);
								}

								//Initialise output vector
								for(i=0;i<ETH_PACKET_LENGTH;i++)
								{
									tmp_float=1.0;
									outvec[i]=*(Xint32*)&tmp_float;
								}
								
								// Send back a dummy UDP packet
								pnew.next = NULL;
								pnew.payload = (unsigned char *)outvec;
								pnew.len = ETH_PACKET_LENGTH*sizeof(Xint32);
								pnew.type = PBUF_RAM;
								pnew.tot_len = pnew.len;
								pnew.ref = 1;
								pnew.flags = 0;

								udp_sendto(pcb, &pnew, addr, port);

					}
					else if (packet_type==3) //write data to SOC input arrays
					{

						if (DEBUG)
						printf("Write data to SOC input arrays ...\n");

						switch (packet_internal_ID)
						{
							case 0: //soc_state0_in
							if (DEBUG)
									printf("write soc_state0_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_STATE0_IN_VECTOR_LENGTH)
										{break;}
									soc_state0_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							case 1: //soc_primal0_in
							if (DEBUG)
									printf("write soc_primal0_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_PRIMAL0_IN_VECTOR_LENGTH)
										{break;}
									soc_primal0_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							case 2: //soc_dual0_in
							if (DEBUG)
									printf("write soc_dual0_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_DUAL0_IN_VECTOR_LENGTH)
										{break;}
									soc_dual0_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							case 3: //soc_tol_iterates_in
							if (DEBUG)
									printf("write soc_tol_iterates_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_TOL_ITERATES_IN_VECTOR_LENGTH)
										{break;}
									soc_tol_iterates_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							default:
							break;
						}
							//Initialise output vector
							for(i=0;i<ETH_PACKET_LENGTH;i++)
							{
								tmp_float=1.0;
								outvec[i]=*(Xint32*)&tmp_float;
							}
							
							// Send back a dummy UDP packet
							pnew.next = NULL;
							pnew.payload = (unsigned char *)outvec;
							pnew.len = ETH_PACKET_LENGTH*sizeof(Xint32);
							pnew.type = PBUF_RAM;
							pnew.tot_len = pnew.len;
							pnew.ref = 1;
							pnew.flags = 0;

							udp_sendto(pcb, &pnew, addr, port);

					}

					else if (packet_type==4) //read data from SOC output vectors
					{
						
						tmp_float=((float)CntValue1-(float)CntValue2) / (float)EE_TICKS_PER_SEC;
						if (DEBUG)
							printf("IP time = %f [s]\n",tmp_float);

						outvec[ETH_PACKET_LENGTH_RECV-1]=*(Xint32*)&tmp_float;  //time


						//Initialise output vector
						for(i=0;i<ETH_PACKET_LENGTH_RECV-2;i++)
						{
							tmp_float=1.0;
							outvec[i]=*(Xint32*)&tmp_float;
						}
						
						
						switch (packet_internal_ID)
						{
							case 0: //primal_in
							if (DEBUG)
								printf("read primal_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_PRIMAL_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_primal_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							case 1: //dual_in
							if (DEBUG)
								printf("read dual_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_DUAL_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_dual_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							case 2: //aux_primal_in
							if (DEBUG)
								printf("read aux_primal_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_AUX_PRIMAL_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_aux_primal_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							case 3: //iterates_in
							if (DEBUG)
								printf("read iterates_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_ITERATES_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_iterates_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							default:
							break;
						}
						
						

						// send back the payload
						// We now need to return the result
						pnew.next = NULL;
						pnew.payload = (unsigned char *)outvec;
						pnew.len = ETH_PACKET_LENGTH*sizeof(Xint32);
						pnew.type = PBUF_RAM;
						pnew.tot_len = pnew.len;
						pnew.ref = 1;
						pnew.flags = 0;


						udp_sendto(pcb, &pnew, addr, port);

					}
					else if (packet_type==5) //read data from DDR
					{

						//Initialise output vector
						for(i=0;i<ETH_PACKET_LENGTH;i++)
						{
							tmp_float=1.0;
							outvec[i]=*(Xint32*)&tmp_float;
						}
						
						// send back the payload
						// We now need to return the result
						pnew.next = NULL;
						pnew.payload = (unsigned char *)outvec;
						pnew.len = ETH_PACKET_LENGTH*sizeof(Xint32);
						pnew.type = PBUF_RAM;
						pnew.tot_len = pnew.len;
						pnew.ref = 1;
						pnew.flags = 0;

						udp_sendto(pcb, &pnew, addr, port);

					}

				}
				/* free the received pbuf */
				pbuf_free(p);

				return ERR_OK;

	//printf("Time was: %d\n", tval);
}


err_t tcp_server_function(void *arg, struct tcp_pcb *tpcb,
                               struct pbuf *p, err_t err)

{

	int k1;
	int i;
	int len;

	int packet_type;
	int packet_internal_ID;
	int packet_internal_ID_offset;

	int timer=0;

	float tmp_float;


	XScuTimer_Config *ConfigPtr;
	XScuTimer *TimerInstancePtr = &Timer;


	xcore.Bus_a_BaseAddress = 0x43c00000;
	xcore.IsReady = XIL_COMPONENT_IS_READY;



	// indicate that the packet has been received
	tcp_recved(tpcb, p->len);

	// Only respond when the packet is the correct length
	if (p->len == (ETH_PACKET_LENGTH)*sizeof(int32_t)){


		/* Pick up a pointer to the payload */
		payload_ptr = (unsigned char *)p->payload;

			//Free the packet buffer
			pbuf_free(p);

			// Get the payload out
			for(k1=0;k1<ETH_PACKET_LENGTH;k1++){
				get_payload();
				inputvec[k1] = payload_temp;
			}

				//extract informations form the input packet
				tmp_float=*(float*)&inputvec[ETH_PACKET_LENGTH-2];
				packet_type=(int)tmp_float & 0x0000FFFF;
				packet_internal_ID=((int)tmp_float & 0xFFFF0000) >> 16; //if write packet_type, packet_num is the data vector ID

				tmp_float=*(float*)&inputvec[ETH_PACKET_LENGTH-1];
				packet_internal_ID_offset=(int)tmp_float; //if write packet_type, packet_num is the data vector ID

			if (DEBUG){
				printf("\n");
				printf("Received packet:\n");
				printf("packet_type=%x\n",packet_type);
				printf("packet_internal_ID=%x\n",packet_internal_ID);
				printf("packet_internal_ID_offset=%x\n",packet_internal_ID_offset);
			}


			if (packet_type==1) //reset IP
						{

							if (DEBUG)
								printf("Reset IP ...\n");

							if (XFoo_IsIdle(&xcore)==1)
							{
								if (DEBUG){
									printf("The core is ready to be used\n");
								}
							}
							else
								printf("ERROR: reprogram the FPGA !\n"); //should be added the IP reset procedure

							//Initialise output vector
							for(i=0;i<ETH_PACKET_LENGTH;i++)
							{
								tmp_float=1.0;
								outvec[i]=*(Xint32*)&tmp_float;
							}

							// send back a dummy packet
							err = tcp_write(tpcb, outvec, ETH_PACKET_LENGTH*sizeof(Xint32), TCP_WRITE_FLAG_MORE);
							tcp_output(tpcb); //send data now
						}
					else if (packet_type==2) //start SOC function
							{

								if (DEBUG)
									printf("Start SOC function ...\n");


								//Start timer
								CntValue1 = IpTimer(TIMER_DEVICE_ID,0);

								//Start SOC function
								soc_user(soc_state0_in,soc_primal0_in,soc_dual0_in,soc_tol_iterates_in,soc_primal_out,soc_dual_out,soc_aux_primal_out,soc_iterates_out);

								//Stop timer
								CntValue2 = IpTimer(TIMER_DEVICE_ID,1);

								if (DEBUG)
								{
									printf ("SOC Timer: beginning of the counter is : %d \n", CntValue1);
									printf ("SOC Timer: end of the counter is : %d \n", CntValue2);
								}
							//Initialise output vector
							for(i=0;i<ETH_PACKET_LENGTH;i++)
							{
								tmp_float=1.0;
								outvec[i]=*(Xint32*)&tmp_float;
							}

							// send back a dummy packet
							err = tcp_write(tpcb, outvec, ETH_PACKET_LENGTH*sizeof(Xint32), TCP_WRITE_FLAG_MORE);
							tcp_output(tpcb); //send data now

					}
					else if (packet_type==3) //write data to SOC input arrays
					{

						if (DEBUG)
						printf("Write data to SOC input arrays ...\n");

						switch (packet_internal_ID)
						{
							case 0: //soc_state0_in
							if (DEBUG)
									printf("write soc_state0_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_STATE0_IN_VECTOR_LENGTH)
										{break;}
									soc_state0_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							case 1: //soc_primal0_in
							if (DEBUG)
									printf("write soc_primal0_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_PRIMAL0_IN_VECTOR_LENGTH)
										{break;}
									soc_primal0_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							case 2: //soc_dual0_in
							if (DEBUG)
									printf("write soc_dual0_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_DUAL0_IN_VECTOR_LENGTH)
										{break;}
									soc_dual0_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							case 3: //soc_tol_iterates_in
							if (DEBUG)
									printf("write soc_tol_iterates_in\n\r");

								for (i=0; i<ETH_PACKET_LENGTH-2; i++)
								{
									tmp_float=*(float*)&inputvec[i];
									if(i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset > SOC_TOL_ITERATES_IN_VECTOR_LENGTH)
										{break;}
									soc_tol_iterates_in[i+(ETH_PACKET_LENGTH-2)*packet_internal_ID_offset] = tmp_float;
								}
							break;
							
							default:
							break;
						}
						//Initialise output vector
						for(i=0;i<ETH_PACKET_LENGTH;i++)
						{
							tmp_float=1.0;
							outvec[i]=*(Xint32*)&tmp_float;
						}

						// send back a dummy packet
						err = tcp_write(tpcb, outvec, ETH_PACKET_LENGTH*sizeof(Xint32), TCP_WRITE_FLAG_MORE);
						tcp_output(tpcb); //send data now

					}

					else if (packet_type==4) //read data from SOC output vectors
					{
						
						tmp_float=((float)CntValue1-(float)CntValue2) / (float)EE_TICKS_PER_SEC;
						if (DEBUG)
							printf("IP time = %f [s]\n",tmp_float);

						outvec[ETH_PACKET_LENGTH_RECV-1]=*(Xint32*)&tmp_float;  //time


						//Initialise output vector
						for(i=0;i<ETH_PACKET_LENGTH_RECV-2;i++)
						{
							tmp_float=1.0;
							outvec[i]=*(Xint32*)&tmp_float;
						}
						
						
						switch (packet_internal_ID)
						{
							case 0: //primal_in
							if (DEBUG)
								printf("read primal_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_PRIMAL_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_primal_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							case 1: //dual_in
							if (DEBUG)
								printf("read dual_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_DUAL_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_dual_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							case 2: //aux_primal_in
							if (DEBUG)
								printf("read aux_primal_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_AUX_PRIMAL_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_aux_primal_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							case 3: //iterates_in
							if (DEBUG)
								printf("read iterates_out\n\r");
							for (i=0; i<ETH_PACKET_LENGTH_RECV-2; i++)
							{
								if (i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset > SOC_ITERATES_OUT_VECTOR_LENGTH)
									break;
								outvec[i]=*(Xint32*)&soc_iterates_out[i+(ETH_PACKET_LENGTH_RECV-2)*packet_internal_ID_offset];
							}
							break;
							
							default:
							break;
						}
						
						

						// send back the payload
						err = tcp_write(tpcb, outvec, ETH_PACKET_LENGTH*sizeof(Xint32), TCP_WRITE_FLAG_MORE);
						tcp_output(tpcb); //send data now

		}
					else if (packet_type==5) //read data from DDR
					{

						//Initialise output vector
						for(i=0;i<ETH_PACKET_LENGTH;i++)
						{
							tmp_float=1.0;
							outvec[i]=*(Xint32*)&tmp_float;
						}

						// send back the payload
						err = tcp_write(tpcb, outvec, ETH_PACKET_LENGTH*sizeof(Xint32), TCP_WRITE_FLAG_MORE);
						tcp_output(tpcb); //send data now

		}

	}

	/* free the received pbuf */
	pbuf_free(p);

	return ERR_OK;
}


err_t tcp_accept_callback(void *arg, struct tcp_pcb *newpcb, err_t err)
{
	static int connection = 1;


		/* set the receive callback for this connection */
		tcp_recv(newpcb, tcp_server_function);


	/* increment for subsequent accepted connections */
	connection++;

	return ERR_OK;
}


int start_application()
{

	err_t err;
	unsigned port = 2007;


	// bind to specified @port
if (TYPE_ETH==1){ //TCP interface
		struct tcp_pcb *pcb;
		pcb = tcp_new();


		//create new TCP PCB structure */
		if (!pcb) {
			xil_printf("Error creating PCB. Out of Memory\n\r");
			return -1;
		}

		err = tcp_bind(pcb, IP_ADDR_ANY, port);

		if (err != ERR_OK) {
			xil_printf("Unable to bind to port %d: err = %d", port, err);
			return -2;
		}

		/* we do not need any arguments to callback functions */
		//tcp_arg(pcb, NULL);

		/* listen for connections */
		pcb = tcp_listen(pcb);
		if (!pcb) {
			xil_printf("Out of memory while tcp_listen\n\r");
			return -3;
		}

		//xil_printf("Got here: start_application\n\r");

		/* specify callback to use for incoming connections */
		tcp_accept(pcb, tcp_accept_callback);

		xil_printf("TCP/IP  FPGA interface framework server started @ port %d\n\r", port);

}

else //UDP interface
{

		struct udp_pcb *pcb;
		pcb = udp_new();

		//create new TCP PCB structure
		if (!pcb) {
			xil_printf("Error creating PCB. Out of Memory\n\r");
			return -1;
		}

		err = udp_bind(pcb, IP_ADDR_ANY, port);

		if (err != ERR_OK) {
			xil_printf("Unable to bind to port %d: err = %d", port, err);
			return -2;
		}

		// specify callback to use for incoming connections
		udp_recv(pcb, udp_server_function, NULL);

		xil_printf("UDP/IP FPGA interface framework server started @ port %d\n\r", port);

}

	return 0;
}
