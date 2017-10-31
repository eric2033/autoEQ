//Portaudio code by M. Farbood
//FFTW base implementation by O. Nieto

#include <stdlib.h>
#include <stdio.h>
#include <portaudio.h>
#include <sndfile.h>
#include <string.h>
#include <ncurses.h>
#include "inputlib.h"
#include <math.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define FRAMES_PER_BUFFER 4096
#define MONO 1
#define STEREO 2
#define SENSITIVITY_INCREMENT 1
#define ATTENUATION_INCREMENT 0.5
#define INITIAL_SENSITIVITY 5
#define INITIAL_ATTENUATION 3

//data struct
typedef struct
{
	float sampleRate;
	float sensitivity;
	float attenuation;
	SNDFILE *infileprimary;
	SF_INFO infileprimaryinfo;
	SNDFILE *infilesecondary;
	SF_INFO infilesecondaryinfo;
} paData;

//global variables
unsigned flags = FFTW_MEASURE;

//wisdom file
const char *filename = "Surgeq.wis";

//global readcount
int readcount = 0;

//initialize sample offset counter for hop size rewinding
int framesleft = FRAMES_PER_BUFFER;
int hop = FRAMES_PER_BUFFER/2;

//check for the start of the signal to read HRTFs and zero OLA buffer
bool start = TRUE;

//window functions
float win[FRAMES_PER_BUFFER];
float coswin[FRAMES_PER_BUFFER/2];

//global pointers for prev arrays for OLA
fftw_complex *primaryprevL = NULL;
fftw_complex *primaryprevR = NULL;
fftw_complex *secondaryprevL = NULL;
fftw_complex *secondaryprevR = NULL;

//callback function declaration
static int paCallback(const void *inputBuffer, void *outputBuffer,
				unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* timeInfo,
				PaStreamCallbackFlags statusFlags,
				void *userData);

int main(int argc, char **argv)
{
	paData data;
	PaStream *stream;
	PaError err;
	PaStreamParameters outputParams;

	//set windows to zero
	memset(win,0,sizeof(float)*FRAMES_PER_BUFFER);
	memset(coswin,0,sizeof(float)*FRAMES_PER_BUFFER/2);

	//main window
	for (int i = 0; i < FRAMES_PER_BUFFER; i++)
	{
		win[i] = sin((M_PI/(float)FRAMES_PER_BUFFER)*(float)i);
	}

	//half cosine window
	for (int i = FRAMES_PER_BUFFER/2; i < FRAMES_PER_BUFFER; i++)
	{
		coswin[i] = win[i];
	}

	//check for input arguments
	if (argc != 3)
	{
		printf ("\nPlease input two audio files as input arguments: \n");
		printf ("	Usage : AutoEQ <filename1> <filename2>\n");
		exit(1);
	}

	//open files
	if ((data.infileprimary = sf_open(argv[1], SFM_READ, &data.infileprimaryinfo)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[1]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	//open files
	if ((data.infilesecondary = sf_open(argv[2], SFM_READ, &data.infilesecondaryinfo)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[2]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	//print file info
	printf("Primary signal audio file: Frames: %d Channels: %d Sample Rate: %d\n",
		(int)data.infileprimaryinfo.frames, data.infileprimaryinfo.channels, data.infileprimaryinfo.samplerate);

	printf("Seconary signal audio file: Frames: %d Channels: %d Sample Rate: %d\n",
		(int)data.infilesecondaryinfo.frames, data.infilesecondaryinfo.channels, data.infilesecondaryinfo.samplerate);

	//check for sample rate
	if (data.infileprimaryinfo.samplerate != data.infilesecondaryinfo.samplerate)
	{
		printf("Error: audio files must have the same sample rate.\nExiting.\n");
		exit(1);
	}

	//init struct parameters
	data.sampleRate = data.infileprimaryinfo.samplerate;
	data.sensitivity = INITIAL_SENSITIVITY;
	data.attenuation = INITIAL_ATTENUATION;

	//initialize and start PortAudio
	err = Pa_Initialize();
	if (err != paNoError)
	{
		printf("PortAudio error: %s\n", Pa_GetErrorText(err));
		printf("\nExiting.\n");
		exit(1);
	}

	//output stream parameters
	outputParams.device = Pa_GetDefaultOutputDevice();
	outputParams.channelCount = STEREO;
	outputParams.sampleFormat = paFloat32;
	outputParams.suggestedLatency =
		Pa_GetDeviceInfo(outputParams.device)->defaultLowOutputLatency;
	outputParams.hostApiSpecificStreamInfo = NULL;

	err = Pa_OpenStream(&stream,
		NULL,
		&outputParams,
		data.sampleRate,
		FRAMES_PER_BUFFER/2,
		paNoFlag,
		paCallback,
		&data);

	if (err != paNoError)
	{
		printf("PortAudio error: open stream: %s\n", Pa_GetErrorText(err));
		exit(2);
	}

	//start audio stream
	err = Pa_StartStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: start stream: %s\n", Pa_GetErrorText(err));
		exit(3);
	}

	//interactive character input
	initscr();
	cbreak();
	noecho();

	char ch;
	ch = '\0';

	mvprintw(0, 0, "Sensitivity: %d dB  Attenuation: %d dB \n \n[a/z] increases/decreases attenuation (dB)\n"\
		"[k/m] increases/decreases sensitivity (dB) \n"\
		"[q] to quit \n", data.sensitivity, data.attenuation);

	while (ch != 'q')
	{
		ch = getch();
		switch(ch){
			case 'k':
			{
				data.sensitivity += SENSITIVITY_INCREMENT;
				break;
			}
			case 'm':
			{
				data.sensitivity -= SENSITIVITY_INCREMENT;
				if (data.sensitivity < 0){
					data.sensitivity = 0;
				}
				break;
			}
			case 'a':
			{
				data.attenuation += ATTENUATION_INCREMENT;
				break;
			}
			case 'z':
			{
				data.attenuation -= ATTENUATION_INCREMENT;
				if (data.attenuation < 0){
					data.attenuation = 0;
				}
				break;
			}
		}
		mvprintw(0, 0, "Sensitivity: %d dB  Attenuation: %d dB \n ", data.sensitivity, data.attenuation);
	}

	//end curses
	endwin();

	//close files
	sf_close(data.infileprimary);
	sf_close(data.infilesecondary);

	//stop stream
	err = Pa_StopStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: stop stream: %s\n", Pa_GetErrorText(err));
	}

	//close stream
	err = Pa_CloseStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: close stream: %s\n", Pa_GetErrorText(err));
	}

	//terminate PortAudio
	err = Pa_Terminate();
	if (err != paNoError)
	{
		printf("PortAudio error: terminate: %s\n", Pa_GetErrorText(err));
	}

	return 0;
}

//callback function
static int paCallback( const void *inputBuffer, void *outputBuffer,
				unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* timeInfo,
				PaStreamCallbackFlags statusFlags,
				void *userData)
{
	//Cast input params
	paData *data = (paData*)userData;
	float *outbuff = (float *)outputBuffer;
	int i;

	//malloc float arrays for reading
	float *realprimarycurr = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER*STEREO);
	float *realsecondarycurr = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER*STEREO);

	//starting 1/2 block
	float *instartprimary = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);
	float *instartsecondary = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);
	fftw_complex *cstartprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *cstartprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *cstartsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *cstartsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);

	//complex input and output buffers
	fftw_complex *cbuffprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *cbuffprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *cbuffsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *cbuffsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);

	//complex arrays for saving
	fftw_complex *savecbuffprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *savecbuffprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *savecbuffsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *savecbuffsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);

	//reset readcount
 	readcount = 0;

	//read audio into buffers
	if (start == TRUE)
	{
		sf_readf_float(data->infileprimary, instartprimary, FRAMES_PER_BUFFER/2);
		sf_readf_float(data->infilesecondary, instartsecondary, FRAMES_PER_BUFFER/2);
		sf_seek(data->infileprimary, 0, SEEK_SET);
		sf_seek(data->infilesecondary, 0, SEEK_SET);
		readcount = sf_readf_float(data->infileprimary, realprimarycurr, FRAMES_PER_BUFFER);
		readcount = sf_readf_float(data->infilesecondary, realsecondarycurr, FRAMES_PER_BUFFER);
	}
	else if (framesleft < hop && start == FALSE)
	{
		sf_seek(data->infileprimary, -(hop-framesleft), SEEK_END);
		sf_seek(data->infilesecondary, -(hop-framesleft), SEEK_END);
		sf_readf_float(data->infileprimary, realprimarycurr, hop-framesleft);
		sf_readf_float(data->infilesecondary, realsecondarycurr, hop-framesleft);
		sf_seek(data->infileprimary, 0, SEEK_SET);
		sf_seek(data->infilesecondary, 0, SEEK_SET);
		sf_readf_float(data->infileprimary, realprimarycurr+(hop-framesleft), FRAMES_PER_BUFFER-(hop-framesleft));
		sf_readf_float(data->infilesecondary, realsecondarycurr+(hop-framesleft), FRAMES_PER_BUFFER-(hop-framesleft));
	}
	else
	{
		sf_seek(data->infileprimary, -hop, SEEK_CUR);
		sf_seek(data->infilesecondary, -hop, SEEK_CUR);
		readcount = sf_readf_float(data->infileprimary, realprimarycurr, FRAMES_PER_BUFFER);
		readcount = sf_readf_float(data->infilesecondary, realsecondarycurr, FRAMES_PER_BUFFER);
	}

	//if signal reached the end, rewind
	if (readcount < FRAMES_PER_BUFFER)
	{
		sf_seek(data->infileprimary, 0, SEEK_SET);
		sf_seek(data->infilesecondary, 0, SEEK_SET);
		framesleft = sf_readf_float(data->infileprimary, realprimarycurr+readcount, FRAMES_PER_BUFFER-readcount);
		framesleft = sf_readf_float(data->infilesecondary, realsecondarycurr+readcount, FRAMES_PER_BUFFER-readcount);
	}

	//fftw plans
	fftw_plan planprimarycurrL, planprimarycurrR, plansecondarycurrL, plansecondarycurrR, planstartsecondaryL, planstartsecondaryR, planstartprimaryL, planstartprimaryR, outplanprimarycurrL, outplanprimarycurrR, outplansecondarycurrL, outplansecondarycurrR, outplanstartprimaryL, outplanstartprimaryR, outplanstartsecondaryL, outplanstartsecondaryR;

	//deinterleave, window, and copy to complex buffers
	if (start == TRUE)
	{
		for (i = 0; i < framesPerBuffer; i++)
		{
			cstartprimaryL[i][0] = coswin[i] * instartprimary[i*2];
			cstartprimaryL[i][1] = 0.f;
			cstartprimaryR[i][0] = coswin[i] * instartprimary[i*2+1];
			cstartprimaryR[i][1] = 0.f;

			cstartsecondaryL[i][0] = coswin[i] * instartsecondary[i*2];
			cstartsecondaryL[i][1] = 0.f;
			cstartsecondaryR[i][0] = coswin[i] * instartsecondary[i*2+1];
			cstartsecondaryR[i][1] = 0.f;
		}
		for (i = 0; i < FRAMES_PER_BUFFER; i++)
		{
			cbuffprimaryL[i][0] =  win[i] * realprimarycurr[i*2];
			cbuffprimaryL[i][1] = 0.f;
			cbuffprimaryR[i][0] =  win[i] * realprimarycurr[i*2+1];
			cbuffprimaryR[i][1] = 0.f;

			cbuffsecondaryL[i][0] =  win[i] * realsecondarycurr[i*2];
			cbuffsecondaryL[i][1] = 0.f;
			cbuffsecondaryR[i][0] =  win[i] * realsecondarycurr[i*2+1];
			cbuffsecondaryR[i][1] = 0.f;
		}
	}
	else
	{
		for (i = 0; i < FRAMES_PER_BUFFER; i++)
		{
			cbuffprimaryL[i][0] =  win[i] * realprimarycurr[i*2];
			cbuffprimaryL[i][1] = 0.f;
			cbuffprimaryR[i][0] =  win[i] * realprimarycurr[i*2+1];
			cbuffprimaryR[i][1] = 0.f;

			cbuffsecondaryL[i][0] =  win[i] * realsecondarycurr[i*2];
			cbuffsecondaryL[i][1] = 0.f;
			cbuffsecondaryR[i][0] =  win[i] * realsecondarycurr[i*2+1];
			cbuffsecondaryR[i][1] = 0.f;
		}
	}

	// //import wisdom file
	// if (fftw_import_wisdom_from_filename(filename) == 0)
	// {
	// 	printf("Error importing wisdom file.");
	// 	exit(1);
	// }

	//create forward plans
	if (start == TRUE)
	{
		planprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryL, outcbuffprimaryL, FFTW_FORWARD, flags);
		planprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryR, outcbuffprimaryR, FFTW_FORWARD, flags);
		plansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryL, outcbuffsecondaryL, FFTW_FORWARD, flags);
		plansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryR, outcbuffsecondaryR, FFTW_FORWARD, flags);
		planstartprimaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartprimaryL, outcstartprimaryL, FFTW_FORWARD, flags);
		planstartprimaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartprimaryR, outcstartprimaryR, FFTW_FORWARD, flags);
		planstartsecondaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartsecondaryL, outcstartsecondaryL, FFTW_FORWARD, flags);
		planstartsecondaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartsecondaryR, outcstartsecondaryR, FFTW_FORWARD, flags);
	}
	else
	{
		planprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryL, outcbuffprimaryL, FFTW_FORWARD, flags);
		planprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryR, outcbuffprimaryR, FFTW_FORWARD, flags);
		plansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryL, outcbuffsecondaryL, FFTW_FORWARD, flags);
		plansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryR, outcbuffsecondaryR, FFTW_FORWARD, flags);
	}

	//export wisdom file
	if (fftw_export_wisdom_to_filename(filename) == 0)
	{
		printf("Error exporting wisdom file.");
		exit(1);
	}

	//execute forward plans
	if (start == TRUE)
	{
		fftw_execute(planprimarycurrL);
		fftw_execute(planprimarycurrR);
		fftw_execute(plansecondarycurrL);
		fftw_execute(plansecondarycurrR);
		fftw_execute(planstartprimaryL);
		fftw_execute(planstartprimaryR);
		fftw_execute(planstartsecondaryL);
		fftw_execute(planstartsecondaryR);
	}
	else
	{
		fftw_execute(planprimarycurrL);
		fftw_execute(planprimarycurrR);
		fftw_execute(plansecondarycurrL);
		fftw_execute(plansecondarycurrR);
	}

	//destroy plans
	if (start == TRUE)
	{
		fftw_destroy_plan(planprimarycurrL);
		fftw_destroy_plan(planprimarycurrR);
		fftw_destroy_plan(plansecondarycurrL);
		fftw_destroy_plan(plansecondarycurrR);
		fftw_destroy_plan(planstartprimaryL);
		fftw_destroy_plan(planstartprimaryR);
		fftw_destroy_plan(planstartsecondaryL);
		fftw_destroy_plan(planstartsecondaryR);
	}
	else
	{
		fftw_destroy_plan(planprimarycurrL);
		fftw_destroy_plan(planprimarycurrR);
		fftw_destroy_plan(plansecondarycurrL);
		fftw_destroy_plan(plansecondarycurrR);
	}

	//AUTOEQ ALGORITHM

	//create inverse plans
	if (start == TRUE)
	{
		outplanprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryL, savecbuffprimaryL, FFTW_BACKWARD, flags);
		outplanprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryR, savecbuffprimaryR, FFTW_BACKWARD, flags);
		outplansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryL, savecbuffsecondaryL, FFTW_BACKWARD, flags);
		outplansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryR, savecbuffsecondaryR,  FFTW_BACKWARD, flags);
		outplanstartprimaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartprimaryL, cstartprimaryL, FFTW_BACKWARD, flags);
		outplanstartprimaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartprimaryR, cstartprimaryR, FFTW_BACKWARD, flags);
		outplanstartsecondaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartsecondaryL, cstartsecondaryL, FFTW_BACKWARD, flags);
		outplanstartsecondaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartsecondaryR, cstartsecondaryR, FFTW_BACKWARD, flags);
	}
	else
	{
		outplanprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryL, savecbuffprimaryL, FFTW_BACKWARD, flags);
		outplanprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryR, savecbuffprimaryR, FFTW_BACKWARD, flags);
		outplansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryL, savecbuffsecondaryL, FFTW_BACKWARD, flags);
		outplansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryR, savecbuffsecondaryR,  FFTW_BACKWARD, flags);
	}

	//execute inverse plans
	if (start == TRUE)
	{
		fftw_execute(outplanprimarycurrL);
		fftw_execute(outplanprimarycurrR);
		fftw_execute(outplansecondarycurrL);
		fftw_execute(outplansecondarycurrR);
		fftw_execute(outplanstartprimaryL);
		fftw_execute(outplanstartprimaryR);
		fftw_execute(outplanstartsecondaryL);
		fftw_execute(outplanstartsecondaryR);
	}
	else
	{		
		fftw_execute(outplanprimarycurrL);
		fftw_execute(outplanprimarycurrR);
		fftw_execute(outplansecondarycurrL);
		fftw_execute(outplansecondarycurrR);
	}

	//destroy plans
	if (start == TRUE)
	{
		fftw_destroy_plan(outplanprimarycurrL);
		fftw_destroy_plan(outplanprimarycurrR);
		fftw_destroy_plan(outplansecondarycurrL);
		fftw_destroy_plan(outplansecondarycurrR);
		fftw_destroy_plan(outplanstartprimaryL);
		fftw_destroy_plan(outplanstartprimaryR);
		fftw_destroy_plan(outplanstartsecondaryL);
		fftw_destroy_plan(outplanstartsecondaryR);
	}
	else
	{
		fftw_destroy_plan(outplanprimarycurrL);
		fftw_destroy_plan(outplanprimarycurrR);
		fftw_destroy_plan(outplansecondarycurrL);
		fftw_destroy_plan(outplansecondarycurrR);
	}

	//window and scale and output real output
	if (start == TRUE)
	{

		for (i = 0; i < framesPerBuffer; i++)
		{
			outbuff[i*2] = (coswin[i] * cstartprimaryL[i][0] / (FRAMES_PER_BUFFER/2)) + (coswin[i] * cstartsecondaryL[i][0] / (FRAMES_PER_BUFFER/2)) + (win[i] * savecbuffprimaryL[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryL[i][0] / (float)FRAMES_PER_BUFFER);
			outbuff[i*2+1] = (coswin[i] * cstartprimaryR[i][0] / (FRAMES_PER_BUFFER/2)) + (coswin[i] * cstartsecondaryR[i][0] / (FRAMES_PER_BUFFER/2)) + (win[i] * savecbuffprimaryR[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryR[i][0] / (float)FRAMES_PER_BUFFER);
		}
	}
	else
	{
		for (i = 0; i < framesPerBuffer/2; i++)
		{
			outbuff[i*2] = (coswin[i] * primaryprevL[i+framesPerBuffer][0] / (float)FRAMES_PER_BUFFER) + (coswin[i] * secondaryprevL[i+framesPerBuffer][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffprimaryL[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryL[i][0] / (float)FRAMES_PER_BUFFER);
			outbuff[i*2+1] = (coswin[i] * primaryprevR[i+framesPerBuffer][0] / (float)FRAMES_PER_BUFFER) + (coswin[i] * secondaryprevR[i+framesPerBuffer][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffprimaryR[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryR[i][0] / (float)FRAMES_PER_BUFFER);
		}
	}

	//point to saved buffers
	primaryprevL = savecbuffprimaryL;
	primaryprevR = savecbuffprimaryR;
	secondaryprevL = savecbuffsecondaryL;
	secondaryprevR = savecbuffsecondaryR;

	//free buffers
	if (start == TRUE)
	{
		fftw_free(cbuffprimaryL);
		fftw_free(cbuffprimaryR);
		fftw_free(outcbuffprimaryL);
		fftw_free(outcbuffprimaryR);
		fftw_free(cbuffsecondaryL);
		fftw_free(cbuffsecondaryR);
		fftw_free(outcbuffsecondaryL);
		fftw_free(outcbuffsecondaryR);
		fftw_free(realprimarycurr);
		fftw_free(realsecondarycurr);

		fftw_free(instartprimary);
		fftw_free(instartsecondary);
		fftw_free(cstartprimaryL);
		fftw_free(cstartprimaryR);
		fftw_free(outcstartprimaryL);
		fftw_free(outcstartprimaryR);
		fftw_free(cstartsecondaryL);
		fftw_free(cstartsecondaryR);
		fftw_free(outcstartsecondaryL);
		fftw_free(outcstartsecondaryR);
		start = FALSE;
	}
	else
	{
		fftw_free(cbuffprimaryL);
		fftw_free(cbuffprimaryR);
		fftw_free(outcbuffprimaryL);
		fftw_free(outcbuffprimaryR);
		fftw_free(cbuffsecondaryL);
		fftw_free(cbuffsecondaryR);
		fftw_free(outcbuffsecondaryL);
		fftw_free(outcbuffsecondaryR);
		fftw_free(realprimarycurr);
		fftw_free(realsecondarycurr);
	}

	return paContinue;
}