BITS 32

WIDTH equ 1920
HEIGHT equ 1080

%ifndef DEBUG
%define FULLSCREEN
%define AUDIO_THREAD
%define GLCHECK
%else
%define NO_AUDIO
%macro GLCHECK 0
	call glGetError
	test eax, eax
	jz %%ok
	int 3
%%ok:
%endmacro
%endif

%ifndef NO_AUDIO
%include "4klang.inc"
extern __4klang_render@4
%else
%define SAMPLE_RATE 44100
%define MAX_SAMPLES 44100*120
%define SAMPLES_PER_TICK 44100/16
%endif

GL_TEXTURE_2D EQU 0x0de1
GL_FRAGMENT_SHADER EQU 0x8b30
GL_UNSIGNED_BYTE EQU 0x1401
;GL_FLOAT EQU 0x1406
GL_RGBA EQU 0x1908
GL_LINEAR EQU 0x2601
GL_TEXTURE_MIN_FILTER EQU 0x2801
;GL_RGBA16F EQU 0x881a
;GL_FRAMEBUFFER EQU 0x8d40
;GL_COLOR_ATTACHMENT0 EQU 0x8ce0

%macro WINAPI_FUNC 2
%if 1
	extern __imp__ %+ %1 %+ @ %+ %2
	%define %1 [__imp__ %+ %1 %+ @ %+ %2]
%else
	extern _ %+ %1 %+ @ %+ %2
	%define %1 _ %+ %1 %+ @ %+ %2
%endif
%endmacro

%ifdef FULLSCREEN
WINAPI_FUNC ChangeDisplaySettingsA, 8
%endif
%ifdef AUDIO_THREAD
WINAPI_FUNC CreateThread, 24
%endif
%ifdef DEBUG
WINAPI_FUNC MessageBoxA, 16
%endif
WINAPI_FUNC ChoosePixelFormat, 8
WINAPI_FUNC CreateWindowExA, 48
WINAPI_FUNC ExitProcess, 4
WINAPI_FUNC GetAsyncKeyState, 4
WINAPI_FUNC GetDC, 4
WINAPI_FUNC PeekMessageA, 20
WINAPI_FUNC SetPixelFormat, 12
WINAPI_FUNC ShowCursor, 4
WINAPI_FUNC SwapBuffers, 4
WINAPI_FUNC waveOutGetPosition, 12
WINAPI_FUNC waveOutOpen, 24
WINAPI_FUNC waveOutWrite, 12
WINAPI_FUNC wglCreateContext, 4
WINAPI_FUNC wglGetProcAddress, 4
WINAPI_FUNC wglMakeCurrent, 8
WINAPI_FUNC glGenTextures, 8
WINAPI_FUNC glBindTexture, 8
WINAPI_FUNC glTexImage2D, 36
WINAPI_FUNC glTexParameteri, 12
WINAPI_FUNC glRects, 16
%ifdef DEBUG
WINAPI_FUNC glGetError, 0
%endif
WINAPI_FUNC CreateCompatibleDC, 4
WINAPI_FUNC CreateDIBSection, 24
WINAPI_FUNC SelectObject, 8
WINAPI_FUNC SetTextColor, 8
WINAPI_FUNC SetBkMode, 8
WINAPI_FUNC CreateFontW, 56
WINAPI_FUNC DrawTextA, 20

%macro FNCALL 1-*
	%rep %0-1
		%rotate -1
		push %1
	%endrep
	%rotate -1
	call %1
%endmacro

%macro GL_FUNC 1
section _ %+ %1 data align=1
%1:
%defstr %[%1 %+ __str] %1
	db %1 %+ __str, 0
%endmacro

GL_FUNC glCreateShaderProgramv
GL_FUNC glUseProgram
GL_FUNC glGetUniformLocation
GL_FUNC glUniform1i
GL_FUNC glUniform1f
;GL_FUNC glGenFramebuffers
;GL_FUNC glBindFramebuffer
;GL_FUNC glFramebufferTexture2D
;GL_FUNC glUniform1fv

%ifdef DEBUG
GL_FUNC glGetProgramInfoLog
%endif

%ifdef DEBUG
	WNDCLASS EQU static_
%else
	%define WNDCLASS 0xc018
%endif

%ifdef FULLSCREEN
section _devmode data align=1
devmode:
	times 9 dd 0
	db 0x9c, 0, 0, 0
	db 0, 0, 0x1c, 0
	times 15 dd 0
	DD	020H, WIDTH, HEIGHT
	times 10 dd 0
%endif

section _pfd data align=1
pfd:
%if 0
	DW	028H, 01H
	DD	025H
	DB	00H, 020H, 00H, 00H, 00H, 00H, 00H, 00H, 08H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 020H, 00H, 00H, 00H, 00H
	DD	00H, 00H, 00H
%else
	DW	00H, 00H
	DD	21H ;025H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	DD	00H, 00H, 00H
%endif

section _binfo data align=1
TEXT_WIDTH EQU 512
TEXT_HEIGHT EQU 1024
bitmapinfo:
	dd bitmapinfo_end - bitmapinfo
	dd TEXT_WIDTH
	dd TEXT_HEIGHT
	dw 1
	dw 32,
	dd 0,
	dd 0, 0, 0, 0, 0
bitmapinfo_end:
	dw 0

section _font data align=1
font:
	dw __utf16__('Consolas'), 0

section _text data align=1
text:
	db 'Notice me sempai', 10
	db 'Farbrausch', 10
	db 'Logicoma', 10
	db 'LJ', 10
	db 'Prismbeings', 10
	db 'Alcatraz', 10
	db 'Conspiracy', 10
	db 'Quite', 10
	db 'Mercury', 10
	db 'SandS', 10
	db 'Titan', 10
	db 'Throb', 10
	db 0

section _rect data align=1
rect:
	dd 0, 0, TEXT_WIDTH, TEXT_HEIGHT
section _wavefmt data align=1
wavefmt:
	dw 3 ; wFormatTag = WAVE_FORMAT_IEEE_FLOAT
	dw 2 ; nChannels
	dd SAMPLE_RATE ; nSamplesPerSec
	dd SAMPLE_RATE * 4 * 2; nAvgBytesPerSec
  dw 4 * 2 ; nBlockAlign
  dw 8 * 4 ; wBitsPerSample
  dw 0 ; cbSize

section _wavehdr data align=1
wavehdr:
	dd sound_buffer ; lpData
	dd MAX_SAMPLES * 2 * 4 ; dwBufferLength
	times 2 dd 0 ; unused stuff
	dd 2 ; dwFlags WHDR_PREPARED  =  0x00000002
	times 4 dd 0 ; unused stuff
	wavehdr_size EQU ($ - wavehdr)

section _mmtime bss align=1
mmtime: resb 12

section _waveout bss align=1
waveout: resd 8

section _sndbuf bss align=1
tmp:
sound_buffer: resd MAX_SAMPLES * 2

%ifdef DEBUG
section _infolog bss align=1
infolog: resb 1024
%endif

section _shader data align=1
%if 1
%include "shader.inc"
%else
_shader_frag:
	db 'uniform int t;'
	db 'float t = t/44100.;'
	db 'void main(){gl_FragColor = vec4(sin(t));}'
%endif

section _shdrptr data align=1
src_main:
	dd _shader_frag

section _strings data align=1
%ifdef DEBUG
static_: db "static", 0
%endif

section _code text align=1
_start:
%if 1
	%define ZERO 0
%else
	%define ZERO ecx
	xor ZERO, ZERO
%endif

%ifdef FULLSCREEN
	FNCALL ChangeDisplaySettingsA, devmode, 4
%endif

	FNCALL ShowCursor, ZERO
	FNCALL CreateWindowExA, ZERO, WNDCLASS, ZERO, 0x90000000, ZERO, ZERO, WIDTH, HEIGHT, ZERO, ZERO, ZERO, ZERO
	FNCALL GetDC, eax
	mov ebp, eax ; ebp = HDC
	FNCALL ChoosePixelFormat, ebp, pfd
	FNCALL SetPixelFormat, ebp, eax, pfd
	FNCALL wglCreateContext, ebp
	FNCALL wglMakeCurrent, ebp, eax
	GLCHECK

%ifndef NO_AUDIO
%ifdef AUDIO_THREAD
	FNCALL CreateThread, ZERO, ZERO, __4klang_render@4, sound_buffer, ZERO, ZERO
%else
	FNCALL __4klang_render@4, sound_buffer
%endif
%endif

	FNCALL CreateCompatibleDC, ZERO
	mov ebx, eax
	FNCALL CreateDIBSection, ebx, bitmapinfo, ZERO, waveout, ZERO, ZERO
	mov edx, eax
	FNCALL SelectObject, ebx, edx
	FNCALL SetTextColor, ebx, 0x00ffffff
	FNCALL SetBkMode, ebx, 1 ;TRANSPARENT
	ANTIALIASED_QUALITY EQU 4
	FNCALL CreateFontW, 44, ZERO, ZERO, ZERO, 700, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ANTIALIASED_QUALITY, ZERO, font
	FNCALL SelectObject, ebx, eax
  ;FNCALL DrawTextA, ebx, text, -1, ZERO, ZERO
	DT_WORDBREAK EQU 0x10
	DT_CALCRECT EQU 0x400
  ;FNCALL DrawTextA, ebx, _shader_frag, -1, rect, (DT_CALCRECT | DT_WORDBREAK)
  ;FNCALL DrawTextA, ebx, _shader_frag, -1, rect, DT_WORDBREAK
  ;FNCALL DrawTextA, ebx, text, -1, rect, DT_CALCRECT
  FNCALL DrawTextA, ebx, text, -1, rect, 0

	TEX_TEXT EQU 1
	FNCALL glGenTextures, 1, mmtime
	FNCALL glBindTexture, GL_TEXTURE_2D, TEX_TEXT
	FNCALL glTexImage2D, GL_TEXTURE_2D, 0, GL_RGBA, TEXT_WIDTH, TEXT_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, dword [waveout]
	FNCALL glTexParameteri, GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR

	FNCALL wglGetProcAddress, glCreateShaderProgramv
	FNCALL eax, GL_FRAGMENT_SHADER, 1, src_main
%ifdef DEBUG
	push eax
	push infolog
	push 0
	push 1023
	push eax
	FNCALL wglGetProcAddress, glGetProgramInfoLog
	call eax
	push 0
	push infolog
	push infolog
	push 0
	call MessageBoxA
	pop eax
%endif
	mov esi, eax
	FNCALL wglGetProcAddress, glUseProgram
	FNCALL eax, esi
	GLCHECK

	FNCALL wglGetProcAddress, glGetUniformLocation
	FNCALL eax, _var_TEX, esi
	mov ebx, eax
	FNCALL wglGetProcAddress, glUniform1i
	FNCALL eax, ebx, 1

	; PLAY MUSIC
	FNCALL waveOutOpen, waveout, byte -1, wavefmt, ZERO, ZERO, ZERO
	FNCALL waveOutWrite, dword [waveout], wavehdr, wavehdr_size

	mainloop:
	;mov ebx, esp
	;mov dword [ebx], 4
	FNCALL waveOutGetPosition, dword [waveout], mmtime, 12
	mov ebx, dword [mmtime + 4]
	cmp ebx, MAX_SAMPLES * 8
	jge exit

	push 01bH ;GetAsyncKeyState

	; PeekMessageA
	push 1
	push 0
	push 0
	push 0
	push 0

	; SwapBuffers
	push ebp

	; glRects
	push 1
	push 1
	push byte -1
	push byte -1

	push ebx
	fild dword [esp]
	push SAMPLES_PER_TICK * 8
	fild dword [esp]
	;mov dword [esp], SAMPLES_PER_TICK * 8
	;fidiv dword [esp]
	fdivp
	pop ebx
	fstp dword [esp]
	push _var_T
	push esi
	push glGetUniformLocation

	call wglGetProcAddress
	call eax
	push eax
	push glUniform1f
	call wglGetProcAddress
	call eax

	call glRects

	call SwapBuffers
	call PeekMessageA
	call GetAsyncKeyState
	jz mainloop
exit:
	call ExitProcess
