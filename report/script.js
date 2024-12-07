document.addEventListener("DOMContentLoaded", function () {
    const introPage = document.getElementById("introPage");
    const mainContent = document.getElementById("mainContent");
    let hasScrolledToTop = false; // 用來追蹤是否已經滾動到頂部

    // 滾動或滾輪事件觸發
    window.addEventListener("wheel", function (event) {
        if (event.deltaY > 0 && !hasScrolledToTop) {
            // 隱藏第一頁
            introPage.style.transform = "translateY(-100%)";
            introPage.style.opacity = "0";

            // 顯示主內容並強制滾動到頂部
            setTimeout(function () {
                mainContent.classList.remove("hidden");
                mainContent.classList.add("visible");
                window.scrollTo({ top: 0, behavior: "smooth" }); // 滾動到頁面頂部
                hasScrolledToTop = true; // 標記已經滾動過一次
            }, 500); // 等待第一頁淡出後顯示
        }
    });
});
