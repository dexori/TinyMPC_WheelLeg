#pragma once
#include <cstdint>

// 由底层工程实现：单调微秒计时，以及保存/恢复原中断状态的临界区。
// 不要无条件开中断，也不要把求解过程整体放进关中断区。
std::uint64_t tinympc_time_us();
std::uint32_t tinympc_enter_critical();
void tinympc_exit_critical(std::uint32_t previous_state);

// 若放入 DTCM，由工程定义本宏并配置对应链接段。
#ifndef TINYMPC_TASK_STORAGE
#define TINYMPC_TASK_STORAGE
#endif
