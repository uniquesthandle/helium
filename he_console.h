// he_console.h - вывод кириллицы в консоли Windows: переключает кодовую страницу консоли на UTF-8.
// Без этого в cmd/PowerShell (cp866/cp1251) русский текст виден как «╨н╤В╨░╨╗╨╛╨╜». На других системах – пустая функция.
// Подключать ПОСЛЕ остальных заголовков (windows.h вводит макросы).
#pragma once
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif
inline void he_console_utf8() {
#if defined(_WIN32)
  SetConsoleOutputCP(65001);
  SetConsoleCP(65001);
#endif
}