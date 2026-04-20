document.addEventListener("DOMContentLoaded", () => {
  const langInner = document.querySelector(".md-header__option .md-select__inner");
  if (!langInner) return;
  const option = langInner.closest(".md-header__option");
  if (!option) return;

  const toTargetPath = (pathname) => {
    if (pathname.includes("/en/")) return pathname.replace("/en/", "/");
    if (pathname.includes("/site/")) return pathname.replace("/site/", "/site/en/");
    if (pathname.startsWith("/en")) return pathname.replace(/^\/en/, "/");
    return `/en${pathname.startsWith("/") ? "" : "/"}${pathname}`;
  };

  option.addEventListener(
    "click",
    (event) => {
      event.preventDefault();
      event.stopPropagation();
      const { pathname, search, hash } = window.location;
      const next = toTargetPath(pathname);
      window.location.assign(`${next}${search}${hash}`);
    },
    true
  );
});
