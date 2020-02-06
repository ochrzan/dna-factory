from multiprocessing import Condition, Value


class SynchCondition:

    def __init__(self, target_thread_count):
        self.target_threads = target_thread_count
        self.current_threads = Value('i', 0)
        self.inner_condition = Condition()

    def wait_for_all(self):
        self.inner_condition.acquire()
        if self.current_threads.value + 1 >= self.target_threads:
            # Release all
            self.current_threads.value = 0
            self.inner_condition.notify_all()
            self.inner_condition.release()
        else:
            self.current_threads.value += 1
            self.inner_condition.wait()
